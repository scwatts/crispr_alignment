#!/usr/bin/env python3
import argparse
import csv
import pathlib


import igraph


class OrderedSpacers():

    def __init__(self, name, spacers):
        self.name = name
        self.spacers = spacers
        self.misorders = list()

    def __str__(self):
        spacers_str = ' '.join(self.spacers)
        misorders_gen = ('({0}, {1})'.format(*e['nodes']) for e in self.misorders)
        misorders_str = ' '.join(misorders_gen) if self.misorders else '-'
        return '\t'.join([self.name, spacers_str, misorders_str])


def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_gffs', nargs='+', required=True, type=pathlib.Path,
            help='Input file path containing multiple CRISPRDetect gffs')
    parser.add_argument('--output_prefix', type=str, default='./crispr_alignment',
            help='Output file prefix [Default: ./crispr_alignment]')

    args = parser.parse_args()

    return args


def main():
    # Get command line arguments
    args = get_arguments()

    # Read in CRISPR data from gff files
    all_spacers = parse_gff_files(args.input_gffs)

    # Create graph and plot
    graph = generate_graph(all_spacers)
    igraph.plot(graph, '%s_graph_plot.png' % args.output_prefix)

    # TODO: pull out neighbourhood graph
    # TODO: assign spacers to one of the graphs, should be relatively clear

    # Pull out subgraphs, assign spacers and then order
    for i, subgraph_nodes in enumerate(graph.clusters(mode=igraph.WEAK), 1):
        subgraph = graph.induced_subgraph(subgraph_nodes, implementation='create_from_scratch')

        # Get spacers which match the subgraph
        subgraph_node_names = [v['name'] for v in list(subgraph.vs)]
        subspacers = dict()
        for spacer_name, spacers in all_spacers.items():
            if len(set(subgraph_node_names).intersection(spacers)) > 0:
                subspacers[spacer_name] = spacers


        # Get output filepath
        output_fp = pathlib.Path('%s_group_%s.tsv' % (args.output_prefix, i))

        # Order spacers via graph
        order_graph_spacers(subgraph, subspacers, output_fp)


def parse_gff_files(input_gffs):
    # Read in spacer CSV
    cas = []
    for gff in input_gffs:
        with open(gff, newline='') as gff_file:
            ca_id = 0
            gff_reader = csv.reader(gff_file, delimiter='\t')
            for row in gff_reader:
                if row[2] == 'repeat_region':
                    ca_id = ca_id + 1
                if row[2] == 'binding_site':
                    ca = [str(gff)]
                    row_info = row[8].split(';')
                    ca[0] = ca[0] + '_' + str(ca_id)
                    for element in row_info:
                        ca.append(element.split('=')[1])
                    cas.append(ca)

    sp_id = 0
    sp_dict = {}

    all_arrays = []
    current_array = []

    name = cas[0][0]
    current_array.append(cas[0][0])

    for ca in cas:
        current_name = ca[0]
        if current_name!=name:
            all_arrays.append(current_array)
            current_array = []
            current_array.append(ca[0])
            name = current_name
        if ca[4] not in sp_dict:
            sp_dict[ca[4]] = str(sp_id)
            sp_id = sp_id + 1
        current_array.append(sp_dict[ca[4]])

    all_arrays.append(current_array)

    all_spacers = {}
    for x in all_arrays:
        all_spacers[x[0]] = x[1:]

    return all_spacers


def order_graph_spacers(graph, all_spacers, output_fp):
    # Get topological order, remove edges to demote graph to acyclic if required
    order_indices, deleted_edge_list = get_spacer_order(graph)
    node_names = list(graph.vs)
    order_names = [node_names[i]['name'] for i in order_indices]

    # If we have deleted edges, collect them into a dict
    deleted_edges = dict()
    if deleted_edge_list:
        for edge in deleted_edge_list:
            try:
                deleted_edges[edge['name']].append(edge)
            except KeyError:
                deleted_edges[edge['name']] = [edge]

    # Create spacer alignment using order indices
    all_ordered_spacers = order_spacers(all_spacers, order_names, deleted_edges)

    # Print output to stdout
    header = ['spacer_name', 'spacer_alignment', 'misordered']
    with output_fp.open('w') as fh:
        print(*header, sep='\t', file=fh)
        for ordered_spacers in all_ordered_spacers:
            print(ordered_spacers, file=fh)


def generate_graph(all_spacers):
    # Create the graph, adding vertice
    graph = igraph.Graph(directed=True)
    vert_set = set()
    for x in all_spacers:
        for y in all_spacers[x]:
            vert_set.add(y)
    for x in sorted(vert_set):
        graph.add_vertex(x)

    # Add vertices from input
    for name, spacers in all_spacers.items():
        for i in range(len(spacers)-1):
            graph.add_edge(spacers[i], spacers[i+1], name=name)

    return graph


def get_spacer_order(graph):
    # Try to run topo sorting. If this fails because graph is cyclic, remove nodes to restore
    # acyclicness
    # TODO: clean this up?
    order_indices = graph.topological_sorting()
    if len(order_indices) != len(graph.vs()):
        # Must preserve and return deleted edges here
        edges_to_delete = graph.feedback_arc_set()
        deleted_edges = collect_edge_info(graph, edges_to_delete)
        graph.delete_edges(edges_to_delete)
        return graph.topological_sorting(), deleted_edges
    else:
        return order_indices, None


def collect_edge_info(graph, edge_indices):
    graph_edges = list(graph.es)
    graph_nodes = list(graph.vs)
    edges = list()
    for edge_index in edge_indices:
        # Get edge data and append to running list
        edge = graph_edges[edge_index]
        edge_data = dict()
        edge_data['name'] = edge['name']
        edge_data['nodes'] = [graph_nodes[i]['name'] for i in edge.tuple]

        edges.append(edge_data)

    return edges


def order_spacers(all_spacers, order_names, deleted_edges):
    # Order spacers lists with determined order
    all_ordered_spacers = list()
    for name, spacers in all_spacers.items():
        order_dict = {s: '-' for s in order_names}
        for spacer in spacers:
            order_dict[spacer] = spacer

        # Init OrderedSpacers instance and add any misorderings
        ordered_spacers = OrderedSpacers(name, order_dict.values())
        if name in deleted_edges:
            ordered_spacers.misorders = deleted_edges[name]
        all_ordered_spacers.append(ordered_spacers)

    return all_ordered_spacers


if __name__ == '__main__':
    main()
