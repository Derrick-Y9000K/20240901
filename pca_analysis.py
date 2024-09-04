from argparse import ArgumentParser
import os
import pandas as pd
import numpy as np
from math import ceil
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from adjustText import adjust_text

def main():
    parser = ArgumentParser()
    parser.add_argument('-d', '--directory', type=str, required=False, help='the file directory')
    parser.add_argument('-i', '--input', type=str, nargs='+', required=True, help='the input file name')
    parser.add_argument('-s', '--separator', type=str, required=False, default='\t', help='the separator used in the matrix file')
    parser.add_argument('-f', '--fontsize', type=int, required=False, default=6, help='the output file name')
    parser.add_argument('-o', '--output', type=str, required=False, help='the output file name')
    parser.add_argument('-3d', '--three-dimensional', action='store_true', default=False, help='the output file')
    args = parser.parse_args()

    if args.directory is not None:
        input = list(map(lambda file: os.path.join(args.directory, file), args.input))
        output = os.path.join(args.directory, 'result.png') if args.output is None else os.path.join(args.directory, args.output)
    else:
        input = args.input
        output = 'result.png' if args.output is None else args.output

    if len(input) == 1:
        data = pd.read_csv(input[0], sep=args.separator, header=0, index_col=0)
    else:
        data = []
        for infile in input:
            data.append(pd.read_csv(infile, sep=args.separator, header=0, index_col=0))
        # data = pd.concat(data, axis=1, join='inner')
        data = pd.concat(data, axis=1, join='outer')
        data.fillna(0, inplace=True)


    data = data.T
    data = np.log2(data + 1)
    data = (data - data.mean()) / data.std()
    pca = PCA(n_components=3, whiten=True)
    result = pca.fit_transform(data)

    group = data.columns.values

    if args.three_dimensional:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.set_xlabel('PC1 (%.2f %%)' % (pca.explained_variance_ratio_[0] * 100), fontdict={'size': 10, 'color': 'black'})
        ax.set_ylabel('PC2 (%.2f %%)' % (pca.explained_variance_ratio_[1] * 100), fontdict={'size': 10, 'color': 'black'})
        ax.set_zlabel('PC3 (%.2f %%)' % (pca.explained_variance_ratio_[2] * 100), fontdict={'size': 10, 'color': 'black'})

        ax.scatter(result[:, 0], result[:, 1], result[:, 2], s=10)
        texts = [ax.text(value[0], value[1], value[2], group[index], fontsize=args.fontsize, color='r') for index, value in enumerate(result)]

    else:
        fig, ax = plt.subplots()
        ax.set_xlabel('PC1 (%.2f %%)' % (pca.explained_variance_ratio_[0] * 100), fontdict={'size': 10, 'color': 'black'})
        ax.set_ylabel('PC2 (%.2f %%)' % (pca.explained_variance_ratio_[1] * 100), fontdict={'size': 10, 'color': 'black'})

        ax.scatter(result[:, 0], result[:, 1], s=10)
        texts = [ax.text(value[0], value[1], group[index], fontsize=args.fontsize, color='r') for index, value in enumerate(result)]

    adjust_text(texts, ax=ax, only_move={'text':'xy'}, arrowprops=dict(arrowstyle='-', lw=1, color='grey'))
    plt.savefig(output, dpi=1200, bbox_inches='tight', pad_inches=0)
    # plt.show()

if __name__ == '__main__':
    main()
