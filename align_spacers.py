#!/usr/bin/env python3
import argparse
import pathlib

import csv
import igraph


def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_gffs', '--list', nargs='+', required=True, type=pathlib.Path,
            help='Input file path containing multiple CRISPRDetect gffs')

    args = parser.parse_args()

    return args


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


def main():
    # Get command line arguments
    args = get_arguments()

    # Read in spacers
    '''
    with args.input_fp.open('r') as fh:
        line_token_gen = (line.rstrip().split() for line in fh)
        all_spacers = {name: spacers for name, *spacers in line_token_gen}
    '''

    # Read in spacer CSV
    cas = []
    for gff in args.input_gffs:
        with open(gff, newline='') as gff_file:
            gff_reader = csv.reader(gff_file, delimiter='\t')
            for row in gff_reader:
                if row[2] == 'binding_site':
                    ca = [gff]
                    row_info = row[8].split(';')
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
        if ca[1].split('_')[0]=='CRISPR1':
            current_name = ca[0]
            if current_name!=name:
                all_arrays.append(current_array)
                current_array = []
                current_array.append(ca[0])
                name = current_name
            if ca[4] not in sp_dict:
                sp_dict[ca[4]] = sp_id
                sp_id = sp_id + 1
            current_array.append(sp_dict[ca[4]])
   
    all_arrays.append(current_array)
    print(*all_arrays,sep='\n')
    
    all_spacers = {}
    for x in all_arrays:
        all_spacers[x[0]] = x[1:]
    
    print(all_spacers[x[0]])
    # Create graph
    graph = generate_graph(all_spacers)

<<<<<<< HEAD
   # Get topological order, remove edges to demote graph to acyclic if required
    order_indices = get_spacer_order(graph)
=======
    # Get topological order, remove edges to demote graph to acyclic if required
    order_indices, deleted_edge_list = get_spacer_order(graph)
>>>>>>> origin/master
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
    print(*header, sep='\t')
    for ordered_spacers in all_ordered_spacers:
        print(ordered_spacers)


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
