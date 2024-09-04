#! /utility/anaconda/envs/RNAseq/bin/python

from argparse import ArgumentParser
import os
import re
import numpy as np
import pandas as pd
from functools import reduce

def main():
    parser = ArgumentParser()
    parser.add_argument('--input-directory', type=str, required=False, default='./', help='the directory with data to merge; default: ./')
    parser.add_argument('--output-directory', type=str, required=False, default='./', help='the directory of results to output; default: ./')
    args = parser.parse_args()
    
    os.makedirs(args.output_directory, exist_ok=True)
    
    result = []
    for file in sorted(os.listdir(args.input_directory)):
        sample = re.search(r'^(?P<name>.*)_count[.]cntTable$', file)['name']
        result.append(pd.read_csv(os.path.join(args.input_directory, file), sep='\t', header=None, names=['gene_id', sample], index_col=0, skiprows=1))
    result = reduce(lambda x, y: pd.merge(x, y, on='gene_id', how='outer'), result).fillna(0)
    
    result.to_csv(os.path.join(args.output_directory, 'TE_count.csv'), header=True, index=True)
    result.to_csv(os.path.join(args.output_directory, 'TE_count.txt'), sep='\t', header=True, index=True)
    
if __name__ == '__main__':
    main()