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


class OrderedSpacers():

    def __init__(self, name, spacers, misordered=False):
        self.name = name
        self.spacers = spacers
        self.misordered = misordered

    def __str__(self):
        misordered_str = '*' if self.misordered else '-'
        spacers_str = ' '.join(self.spacers)
        return '\t'.join([self.name, spacers_str, misordered_str])


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
    all_ordered_spacers = order_spacers(all_spacers, order_names)

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
        graph.delete_edges(graph.feedback_arc_set())
        return graph.topological_sorting()
    else:
        return order_indices


def order_spacers(all_spacers, order_names):
    # Order spacers lists with determined order
    all_ordered_spacers = list()
    for name, spacers in all_spacers.items():
        order_dict = {s: '-' for s in order_names}
        for spacer in spacers:
            order_dict[spacer] = spacer

        # Check if the spacer retains biological order and init OrderedSpacer object
        if spacers == [n for n in order_dict.values() if n != '-']:
            ordered_spacers = OrderedSpacers(name, order_dict.values())
        else:
            ordered_spacers = OrderedSpacers(name, order_dict.values(), misordered=True)
        all_ordered_spacers.append(ordered_spacers)

    return all_ordered_spacers


if __name__ == '__main__':
    main()
