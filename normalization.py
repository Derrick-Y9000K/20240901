import os
import re
import pandas as pd
import numpy as np
from bioinfokit.analys import norm
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('-d', '--directory', type=str, required=False, help='the file directory')
parser.add_argument('-i', '--input', type=str, required=False, default='gene_count_matrix.csv', help='the input gene count matrix file name; default: gene_count_matrix.csv')
parser.add_argument('-l', '--length', type=str, required=False, default='gene_length.txt', help='the gene length file name; default: gene_length.txt')
parser.add_argument('-s', '--separator', type=str, required=False, default=',', help='the separator used in the matrix file; default: ,')
parser.add_argument('-r', '--rpkm-result', type=str, required=False, default='rpkm_matrix.csv', help='the output rpkm matrix file name; default: rpkm_matrix.csv')
parser.add_argument('-t', '--tpm-result', type=str, required=False, default='tpm_matrix.csv', help='the output tpm matrix file name; default: tpm_matrix.csv')
args = parser.parse_args()

def main():
    if args.directory is not None:
        input = os.path.join(args.directory, args.input)
        rpkm_output = os.path.join(args.directory, args.rpkm_result)
        tpm_output = os.path.join(args.directory, args.tpm_result)
    else:
        input = args.input
        rpkm_output = args.rpkm_result
        tpm_output = args.tpm_result

    data = pd.read_csv(input, sep=args.separator, header=0)
    data = data[data.apply(lambda value: '|' in value['gene_id'], axis=1)]
    data['gene_id'] = list(map(lambda value: re.search('^.*[|](?P<gene_name>.*)$', value)['gene_name'], data['gene_id']))
    data.set_index('gene_id', inplace=True)
    data = data[~data.index.duplicated()]

    length = pd.read_csv(args.length, sep='\t', header=None, index_col=0, names=['symbol', 'gene_id', 'region', 'length'])
    length = length[~length.index.duplicated()]
    length.reset_index(drop=False, inplace=True)
    result = pd.merge(data, length, how='left', left_on='gene_id', right_on='symbol')
    result.set_index('symbol', inplace=True)
    result.drop(['gene_id', 'region'], axis=1, inplace=True)

    normalizer = norm()
    normalizer.rpkm(df=result, gl='length')
    rpkm = normalizer.rpkm_norm
    normalizer.tpm(df=result, gl='length')
    tpm = normalizer.tpm_norm

    rpkm.to_csv(rpkm_output, sep=args.separator)
    tpm.to_csv(tpm_output, sep=args.separator)


if __name__ == '__main__':
    main()
