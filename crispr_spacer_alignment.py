#!/usr/bin/env python3
import argparse
import json
import pathlib
import re
import shutil
import subprocess
import sys
import tempfile


import igraph


RC_TABLE = str.maketrans('atgcATGC', 'tacgTACG')
CDHIT_RE = re.compile('^.+>([0-9]+).+$')


class Crispr:

    def __init__(self, contig, start, end, spacers_seqs):
        self.contig = contig
        self.start = start
        self.end = end
        self.spacers_seqs = spacers_seqs

        self.spacers = list()


class OrderedSpacers:

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
    parser.add_argument('--input_fps', nargs='+', required=True, type=pathlib.Path,
            help='Input file path containing multiple CRISPRDetect filepaths')
    parser.add_argument('--output_prefix', type=str, default='./crispr_alignment',
            help='Output file prefix [Default: ./crispr_alignment]')

    args = parser.parse_args()

    return args


def main():
    # Get command line arguments
    args = get_arguments()

    # Read in CRISPR data from json files
    crisprs = parse_json_files(args.input_fps)

    # Cluster spacer sequences using CD-HIT
    spacers_clusters = cluster_spacer_sequences(crisprs)

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


def parse_json_files(input_fps):
    crisprs = list()
    for input_fp in input_fps:
        with input_fp.open('r') as fh:
            input_data = json.loads(fh.read())
            date, version, command, contigs_data = input_data.values()
            for crispr in collect_crispr_from_json(contigs_data):
                crisprs.append(crispr)
    return crisprs


def collect_crispr_from_json(contigs_data):
    for contig_data in contigs_data:
        # Skip contigs without crisprs
        if not contig_data['Crisprs']:
            continue

        crispr_data_gen = (d for d in contig_data['Crisprs'])
        contig_name = contig_data['Id']

        # Extract each crispr found on this contig
        for crispr_data in crispr_data_gen:
            start = crispr_data['Start']
            end = crispr_data['End']
            spacers = [r['Sequence'] for r in crispr_data['Regions'] if r['Type'] == 'Spacer']
            crispr = Crispr(contig_name, start, end, spacers)
            yield crispr


def cluster_spacer_sequences(crisprs):
    # Get spacer sequences and find canonically smallest for each
    canon_spacers = set()
    spacer_gen = (s for c in crisprs for s in c.spacers_seqs)
    for spacer in spacer_gen:
        spacer_rc = spacer[::-1].translate(RC_TABLE)
        if spacer <= spacer_rc:
            canon_spacers.add(spacer)
        else:
            canon_spacers.add(spacer_rc)

    # Convert canon_spacers to list, downstream processing requires
    canon_spacers = list(canon_spacers)

    with tempfile.TemporaryDirectory() as dh:
        # Write out sequences as FASTA
        canon_spacers_fp = pathlib.Path(dh, 'spacers.fasta')
        with canon_spacers_fp.open('w') as fh:
            for i, spacer in enumerate(canon_spacers, 0):
                fh.write('>%s\n' % i)
                fh.write('%s\n' % spacer)

        # Cluster and parse output
        output_fasta_fp = pathlib.Path(dh, 'result.fasta')
        output_clusters_fp = pathlib.Path(dh, 'result.fasta.clstr')
        run_cdhit(canon_spacers_fp, output_fasta_fp)

        with output_clusters_fp.open('r') as fh:
            cluster_data = fh.readlines()
    spacer_clusters_indices = parse_cdhit_clusters(cluster_data)

    # Create cluster sets containing spaces with rc
    spacers_clusters = list()
    for index_set in  spacer_clusters_indices:
        spacer_set = {canon_spacers[i] for i in index_set}
        spacer_rc_set = {s[::-1].translate(RC_TABLE) for s in spacer_set}
        spacers_clusters.append(spacer_set | spacer_rc_set)
    return spacers_clusters


def run_cdhit(input_fp, output_fp):
    executable = 'cd-hit'

    command = '%s -i %s -o %s -c 0.95 -d 0' % (executable, input_fp, output_fp)

    # Check for CD-HIT in path
    if not shutil.which(executable):
        print('error: could not find %s in PATH' % executable, file=sys.stderr)
        sys.exit(1)

    # Run CD-HIT
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            shell=True, encoding='utf-8')

    # Check retcode before returning
    if result.returncode != 0:
        print('error: failed to run CD-HIT:', file=sys.stderr)
        print(result.stderr, file=sys.stderr)
        sys.exit(1)
    return result


def parse_cdhit_clusters(data):
    clusters = list()
    cluster_members = set()
    for line in data:
        if line.startswith('>'):
            if cluster_members:
                clusters.append(cluster_members)
                cluster_members = set()
            continue
        print(line)
        cluster_number = int(CDHIT_RE.match(line).group(1))
        cluster_members.add(cluster_number)
    return clusters


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
