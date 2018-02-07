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

    # Try to run topo sorting. If this fails because graph is cyclic, remove nodes to restor
    # acyclicness
    # TODO: clean this up
    ordering = graph.topological_sorting()
    if len(ordering) != len(graph.vs()):
        rnodes = graph.feedback_arc_set()
        graph.delete_edges(graph.feedback_arc_set())
    ordering = graph.topological_sorting()

    # Order spacer lists with determined order
    # TODO: see output of ordering and then apply
    print(ordering)
    print(*(list(graph.vs)[i]['name'] for i in ordering))
    print(*graph.get_adjacency(), sep='\n')

    # Print out results, for now


if __name__ == '__main__':
    main()
