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


    # TODO: consider node input space and add prior
    # Generate graph structure
    graph = igraph.Graph(directed=True)
    for name, spacers in all_spacers:
        # Get all edge pairs
        spacer_iter = (spacers[i:i+2] for i in range(len(spacers)))
        edge_gen = (e for e in spacer_iter if len(e) > 1)

        # Add to graph
        for source_node, target_node for edge_gen:
            try:
                graph.add_edge(source_node, target_node)
            except igraph._igraph.InternalError:
                if
                graph.add_vertex(

    # Try to run topo sorting

    # Otherwise approximate feedback arc set

    # Order spacer lists with determined order

    # Print out results, for now


if __name__ == '__main__':
    main()
