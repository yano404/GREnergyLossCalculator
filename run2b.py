#!/bin/env python3

import pandas as pd
import argparse

def run2b(run, run2b_file):
    # Load run2b_file
    df_b = pd.read_csv(
        run2b_file,
        delim_whitespace=True,
        comment='#',
        index_col=0,
        names=['b'])
    # Get magnetic field
    b = df_b.loc[run].b
    return b

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description='Return GR magnetic field')
    parser.add_argument(
            'run',
            type=int,
            help='Run number')
    parser.add_argument(
            '-f',
            '--file',
            default='b.dat',
            help='List of magnetic field')

    args = parser.parse_args()

    # Run number
    run = args.run
    run2b_file = args.file

    print(run2b(run, run2b_file))

