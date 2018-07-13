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

        self.name = '%s_%s_%s' % (self.contig, self.start, self.end)

        self.spacers = list()
        self.strong_order = dict()
        self.strong_misorders = list()


def get_arguments():
    parser_parent = argparse.ArgumentParser()
    parser_io = parser_parent.add_argument_group('Input and output')
    parser_cluster = parser_parent.add_argument_group('Clustering')
    parser_cdhit = parser_parent.add_argument_group('CD-HIT parameters')

    parser_io.add_argument('--input_fps', nargs='+', required=True, type=pathlib.Path,
            help='Input file path containing multiple CRISPRDetect filepaths')
    parser_io.add_argument('--output_prefix', type=str, default='./crispr_alignment',
            help='Output file prefix [Default: ./crispr_alignment]')

    parser_cluster.add_argument('--reverse_complement', action='store_true', default=True,
            help='Consider reverse complement of sequence during clustering [default: True]')
    parser_cluster.add_argument('--absoulte_identity', action='store_true', default=False,
            help='Use absolute identity for clustering, skip CD-HIT [default: False]')

    parser_cdhit.add_argument('--identity', type=str, default='0.90',
            help='Sequence identity threshold [default: 0.90]')

    args = parser_parent.parse_args()

    return args


def main():
    # Require Python 3.6 for set, dict order guarantee
    if sys.version_info < (3, 6):
        print('error: %s requires Python 3.6 or above' %  pathlib.Path(__file__).name)
        sys.exit(1)

    # Get command line arguments
    args = get_arguments()
    print(args)

    # Read in CRISPR data from json files
    crisprs = parse_json_files(args.input_fps)

    # Cluster spacer sequences
    spacers_clusters = cluster_spacer_sequences(crisprs, args)

    # Assign spacers an appropriate cluster identifier
    for crispr in crisprs:
        assign_spacer_clusters(crispr, spacers_clusters)

    # Create graph and plot
    graph = create_spacer_graph(crisprs)
    igraph.plot(graph, '%s_graph_plot.png' % args.output_prefix)

    # Order spacers via graph
    order_graph_spacers(graph, crisprs)

    # Write out strong spacer order
    output_fp = pathlib.Path('%s.tsv' % args.output_prefix)
    header = ['spacer_name', 'start', 'end', 'spacer_alignment', 'misordered']
    with output_fp.open('w') as fh:
        print(*header, sep='\t', file=fh)
        for crispr in crisprs:
            print(crispr.contig, crispr.start, crispr.end, sep='\t', end='\t', file=fh)
            print(*crispr.strong_order.values(), sep=' ', end='\t', file=fh)
            misorders = ('(%s, %s)' % (a, b) for a, b in crispr.strong_misorders)
            print(*misorders, sep=', ', file=fh)


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


def cluster_spacer_sequences(crisprs, args):
    if args.reverse_complement:
        spacers = canonical_spacers(crisprs)
    else:
        spacers = [s for c in crisprs for s in c.spacers_seqs]

    if args.absoulte_identity:
        pass
    else:
        with tempfile.TemporaryDirectory() as dh:
            # Write out sequences as FASTA
            spacers_fp = pathlib.Path(dh, 'spacers.fasta')
            with spacers_fp.open('w') as fh:
                for i, spacer in enumerate(spacers, 0):
                    fh.write('>%s\n' % i)
                    fh.write('%s\n' % spacer)

            # Cluster and parse output
            output_fasta_fp = pathlib.Path(dh, 'result.fasta')
            output_clusters_fp = pathlib.Path(dh, 'result.fasta.clstr')
            word_lendth = cdhit_word_length(float(args.identity))
            run_cdhit(spacers_fp, output_fasta_fp, args.identity, args.word_length)

            with output_clusters_fp.open('r') as fh:
                cluster_data = fh.readlines()

        spacer_clusters_indices = parse_cdhit_clusters(cluster_data)

    # Create cluster sets containing spacers, using str(i) for operability with igraph
    spacers_clusters = dict()
        for i, index_set in  enumerate(spacer_clusters_indices, 1):
            spacer_set = {canon_spacers[j] for j in index_set}
            if args.reverse_complement:
                spacer_rc_set = {s[::-1].translate(RC_TABLE) for s in spacer_set}
                spacers_clusters[str(i)] = spacer_set | spacer_rc_set
            else:
                spacers_clusters[str(i)] = spacer_set
    return spacers_clusters


def canonical_spacer_sequence(crisprs):
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
    return list(canon_spacers)


def run_cdhit(input_fp, output_fp, identity, word_length):
    executable = 'cd-hit'
    command_template = '%s -i %s -o %s -c %s -n %s -d 0'
    command = command_template % (executable, input_fp, output_fp, identity, word_length)

    # Run CD-HIT and check retcode
    if not shutil.which(executable):
        print('error: could not find %s in PATH' % executable, file=sys.stderr)
        sys.exit(1)

    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            shell=True, encoding='utf-8')

    if result.returncode != 0:
        print('error: failed to run CD-HIT:', file=sys.stderr)
        print(result.stderr, file=sys.stderr)
        sys.exit(1)
    return result


def cdhit_word_length(identity):
    if identity >= 0.7:
        return 5
    elif identity >= 0.6 and identity < 0.7:
        return 4
    elif identity >= 0.5 and identity < 0.6:
        return 3
    elif identity >= 0.4 and identity < 0.5:
        return 2


def parse_cdhit_clusters(data):
    clusters = list()
    cluster_members = set()
    for line in data:
        if line.startswith('>'):
            if cluster_members:
                clusters.append(cluster_members)
                cluster_members = set()
            continue
        cluster_number = int(CDHIT_RE.match(line).group(1))
        cluster_members.add(cluster_number)
    clusters.append(cluster_members)
    return clusters


def assign_spacer_clusters(crispr, spacers_clusters):
    # TODO: quantify bottleneck
    for spacer_seq in crispr.spacers_seqs:
        for spacer_set_id, spacer_set in spacers_clusters.items():
            if spacer_seq in spacer_set:
                crispr.spacers.append(spacer_set_id)
                break
        else:
            raise ValueError('could not find appropriate spacer set')


def create_spacer_graph(crisprs):
    # Create the graph, adding vertices and edges
    graph = igraph.Graph(directed=True)
    spacers = {s for c in crisprs for s in c.spacers}
    for vertex in sorted(spacers):
        graph.add_vertex(vertex)

    for crispr in crisprs:
        for i in range(len(crispr.spacers) - 1):
            graph.add_edge(crispr.spacers[i], crispr.spacers[i+1], name=crispr.name)
    return graph


def order_graph_spacers(graph, crisprs):
    # Get topological order, remove edges to demote graph to acyclic if required
    order_indices, deleted_edges = get_spacer_order(graph)
    node_names = list(graph.vs)
    order_names = [node_names[i]['name'] for i in order_indices]

    # Create spacer alignment using order indices
    order_spacers(crisprs, order_names, deleted_edges)


def get_spacer_order(graph):
    # Run topo sorting. If graph is cyclic, remove nodes to restore acyclicness
    order_indices = graph.topological_sorting()
    if len(order_indices) != len(graph.vs()):
        edges_to_delete = graph.feedback_arc_set()
        deleted_edges = collect_edge_info(graph, edges_to_delete)
        graph.delete_edges(edges_to_delete)
        return graph.topological_sorting(), deleted_edges
    else:
        return order_indices, None


def collect_edge_info(graph, edge_indices):
    graph_edges = list(graph.es)
    graph_nodes = list(graph.vs)
    edges = dict()
    for edge_index in edge_indices:
        edge = graph_edges[edge_index]
        contig_name = edge['name']
        nodes = [graph_nodes[i]['name'] for i in edge.tuple]
        try:
            edges[contig_name].append((nodes))
        except KeyError:
            edges[contig_name] = [(nodes)]
    return edges


def order_spacers(crisprs, order_names, deleted_edges):
    # Provide and fill out strong ordering to Crispr instances
    all_ordered_spacers = list()
    for crispr in crisprs:
        crispr.strong_order = {s: '-' for s in order_names}
        for spacer in crispr.spacers:
            crispr.strong_order[spacer] = spacer

        # Record incongruencies between strong and biological ordering
        if deleted_edges and crispr.name in deleted_edges:
            crispr.strong_misorders = deleted_edges[crispr.name]


if __name__ == '__main__':
    main()
