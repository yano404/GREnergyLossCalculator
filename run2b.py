#!/bin/env python3

import pandas as pd
import argparse

df_b = pd.read_csv(
        'b.dat',
        delim_whitespace=True,
        comment='#',
        index_col=0,
        names=['b'])

def run2b(run):
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
    args = parser.parse_args()

    # Run number
    run = args.run

    print(run2b(run))

