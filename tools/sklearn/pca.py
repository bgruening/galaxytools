import argparse
import numpy as np
from sklearn.decomposition import PCA, IncrementalPCA, KernelPCA


def main():
    parser = argparse.ArgumentParser(description='RDKit screen')
    parser.add_argument('-i', '--infile',
                        help="Input file")
    parser.add_argument('-if', '--informat',
                        help="Input format.")
    parser.add_argument('--header', action='store_true', help="Include the header row or skip it")
    parser.add_argument('-c', '--columns', type=str.lower, default='all', choices=['all', 'include', 'exclude'],
                        help="Choose to select all columns, or exclude/include some")
    parser.add_argument('-ci', '--column_indices', type=str.lower,
                        help="Choose to select all columns, or exclude/include some")
    parser.add_argument('-n', '--number', nargs='?', type=int, default=None,\
                        help="Number of components to keep. If not set, all components are kept")
    parser.add_argument('--whiten', action='store_true', help="Whiten the components")
    parser.add_argument('-t', '--pca_type', type=str.lower, default='classical', choices=['classical', 'incremental', 'kernel'],
                        help="Choose which flavour of PCA to use")
    parser.add_argument('-s', '--svd_solver', type=str.lower, default='auto', choices=['auto', 'full', 'arpack', 'randomized'],
                        help="Choose the type of svd solver.")
    parser.add_argument('-b', '--batch_size', nargs='?', type=int, default=None,\
                        help="The number of samples to use for each batch")
    parser.add_argument('-k', '--kernel', type=str.lower, default='linear',\
                        choices=['linear', 'poly', 'rbf', 'sigmoid', 'cosine', 'precomputed'],
                        help="Choose the type of kernel.")
    parser.add_argument('-g', '--gamma', nargs='?', type=float, default=None,
                        help='Kernel coefficient for rbf, poly and sigmoid kernels. Ignored by other kernels')
    parser.add_argument('-tol', '--tolerance', type=float, default=0.0,
                        help='Convergence tolerance for arpack. If 0, optimal value will be chosen by arpack')
    parser.add_argument('-mi', '--max_iter', nargs='?', type=int, default=None,\
                        help="Maximum number of iterations for arpack")
    parser.add_argument('-d', '--degree', type=int, default=3,\
                        help="Degree for poly kernels. Ignored by other kernels")
    parser.add_argument('-cf', '--coef0', type=float, default=1.0,
                        help='Independent term in poly and sigmoid kernels')
    parser.add_argument('-e', '--eigen_solver', type=str.lower, default='auto', choices=['auto', 'dense', 'arpack'],
                        help="Choose the type of eigen solver.")
    parser.add_argument('-o', '--outfile',
                        help="Base name for output file (no extension).")
    parser.add_argument('-of', '--outformat',
                        help="Output format")
    args = parser.parse_args()

    usecols = None
    pca_params = {}

    with open(args.infile) as f:
        num_cols = len(f.readline().split(' '))

    if args.columns != 'all':
        cols = [int(i) for i in args.column_indices.split(' ')]


    if args.columns == "exclude":
        usecols = sorted(set(range(num_cols)) - set(cols))
    elif args.columns == "include":
        usecols = sorted(set(range(num_cols)).intersection(set(cols)))

    pca_input = np.loadtxt(fname=args.infile, skiprows=int(args.header), usecols=usecols)

    pca_params.update({'n_components': args.number})

    if args.pca_type == 'classical':
        pca_params.update({'svd_solver': args.svd_solver, 'whiten': args.whiten})
        if args.svd_solver == 'arpack':
            pca_params.update({'tol': args.tolerance})
        pca = PCA()

    elif args.pca_type == 'incremental':
        pca_params.update({'batch_size': args.batch_size, 'whiten': args.whiten})
        pca = IncrementalPCA()

    elif args.pca_type == 'kernel':
        pca_params.update({'kernel': args.kernel, 'eigen_solver': args.eigen_solver})

        if args.kernel == 'poly':
            pca_params.update({'gamma': args.gamma, 'degree': args.degree, 'coef0': args.coef0})
        elif args.kernel == 'rbf':
            pca_params.update({'gamma': args.gamma})
        elif args.kernel == 'sigmoid':
            pca_params.update({'gamma': args.gamma, 'coef0': args.coef0})

        if args.eigen_solver == 'arpack':
            pca_params.update({'tol': args.tolerance, 'max_iter': args.max_iter})

        pca = KernelPCA()

    pca.set_params(**pca_params)
    pca_output = pca.fit_transform(pca_input)
    np.savetxt(args.outfile, pca_output)


if __name__ == "__main__":
    main()