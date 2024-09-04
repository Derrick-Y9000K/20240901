import argparse
import os
import re
import numpy as np
import pandas as pd

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, required=True, help='the input gene count matrix file to filter')
    parser.add_argument('-o', '--output', type=str, required=True, help='the output filtered gene count matrix file')
    args = parser.parse_args()
    
    data = pd.read_csv(args.input, header=0, index_col=0)
    data = data.loc[['|' in gene_id for gene_id in data.index]]
    data.index = [gene_id.split('|')[1] for gene_id in data.index]
    data.to_csv(args.output, sep='\t')

if __name__ == '__main__':
    main()