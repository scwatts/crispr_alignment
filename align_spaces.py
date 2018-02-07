#!/usr/bin/env python3
import argparse
import pathlib


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
        spacers = {name: spacers for name, *spacers in line_token_gen}

    # Convert spacer lists to adjacency

    # Try to run topo sorting

    # Otherwise approximate feedback arc set

    # Order spacer lists with determined order

    # Print out results, for now


if __name__ == '__main__':
    main()
