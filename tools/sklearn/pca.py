import argparse
import numpy as np
from sklearn.decomposition import PCA, IncrementalPCA, KernelPCA
from galaxy_ml.utils import read_columns

def main():
    parser = argparse.ArgumentParser(description='RDKit screen')
    parser.add_argument('-i', '--infile',
                        help="Input file")
    parser.add_argument('--header', action='store_true', help="Include the header row or skip it")
    parser.add_argument('-c', '--columns', type=str.lower, default='all', choices=['by_index_number', 'all_but_by_index_number',\
                        'by_header_name', 'all_but_by_header_name', 'all_columns'],
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
    args = parser.parse_args()

    usecols = None
    cols = []
    pca_params = {}

    if args.columns == 'by_index_number' or args.columns == 'all_but_by_index_number':
        usecols = [int(i) for i in args.column_indices.split(',')]
    elif args.columns == 'by_header_name' or args.columns == 'all_but_by_header_name':
        usecols = args.column_indices

    pca_input = read_columns(
        f=args.infile,
        c=usecols,
        c_option=args.columns,
        sep='\t',
        header=int(args.header),
        parse_dates=True,
        encoding=None,
        index_col=None)

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
        pca_params.update({'kernel': args.kernel, 'eigen_solver': args.eigen_solver, 'gamma': args.gamma})

        if args.kernel == 'poly':
            pca_params.update({'degree': args.degree, 'coef0': args.coef0})
        elif args.kernel == 'sigmoid':
            pca_params.update({'coef0': args.coef0})
        elif args.kernel == 'precomputed':
            pca_input = np.dot(pca_input, pca_input.T)

        if args.eigen_solver == 'arpack':
            pca_params.update({'tol': args.tolerance, 'max_iter': args.max_iter})

        pca = KernelPCA()

    pca.set_params(**pca_params)
    pca_output = pca.fit_transform(pca_input)
    np.savetxt(fname=args.outfile, X=pca_output, fmt='%.4f', delimiter='\t')


if __name__ == "__main__":
    main()