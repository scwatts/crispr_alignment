#!/usr/bin/env python3
import argparse
import pathlib


import igraph


def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_fp', required=True, type=pathlib.Path,
            help='Input file path containing crispr spacers')

    args = parser.parse_args()
    if not args.input_fp.exists():
        parser.error('Input file %s does not exist' % args.input_fp)

    return args


def main():
    # Get command line arguments
    args = get_arguments()

    # Read in spacers
    with args.input_fp.open('r') as fh:
        line_token_gen = (line.rstrip().split() for line in fh)
        all_spacers = {name: spacers for name, *spacers in line_token_gen}

    # Create graph
    graph = generate_graph(all_spacers)

    # Get topological order, remove edges to demote graph to acyclic if required
    order_indices = get_spacer_order(graph)
    node_names = list(graph.vs)
    order_names = [node_names[i]['name'] for i in order_indices]

    # Create spacer alignment using order indices
    alignment = create_alignment(all_spacers, order_names)
    for name, ordered_spacers in alignment.items():
        print(name, *ordered_spacers, sep='\t')


def generate_graph(all_spacers):
    # Create the graph, adding vertice
    graph = igraph.Graph(directed=True)
    vert_set = set()
    for x in all_spacers:
        for y in all_spacers[x]:
            vert_set.add(y)
    for x in vert_set:
        graph.add_vertex(x)

    # Add vertices from input
    for x in all_spacers:
        for i in range(len(all_spacers[x])-1):
            graph.add_edge(all_spacers[x][i],all_spacers[x][i+1])

    return graph


def get_spacer_order(graph):
    # Try to run topo sorting. If this fails because graph is cyclic, remove nodes to restore
    # acyclicness
    # TODO: clean this up?
    order_indices = graph.topological_sorting()
    if len(order_indices) != len(graph.vs()):
        rnodes = graph.feedback_arc_set()
        graph.delete_edges(graph.feedback_arc_set())
        return graph.topological_sorting()
    else:
        return order_indices


def create_alignment(all_spacers, order_names):
    # Order spacers lists with determined order
    ordered_spacers = dict()
    for name, spacers in all_spacers.items():
        order_dict = {s: '-' for s in order_names}
        for spacer in spacers:
            order_dict[spacer] = spacer
        ordered_spacers[name] = list(order_dict.values())
    return ordered_spacers


if __name__ == '__main__':
    main()