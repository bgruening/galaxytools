"""Dataset handler."""

import numpy as np

import pandas as pd

__author__ = "Gianluca Corrado"
__copyright__ = "Copyright 2016, Gianluca Corrado"
__license__ = "MIT"
__maintainer__ = "Gianluca Corrado"
__email__ = "gianluca.corrado@unitn.it"
__status__ = "Production"


class Dataset(object):
    """General dataset."""

    def __init__(self, fp, fr, standardize_proteins=False,
                 standardize_rnas=False):
        """
        Constructor.

        Parameters
        ----------
        fp : str
            Protein features

        fr : str
            The name of the HDF5 file containing features for the RNAs.
        """
        self.Fp = fp.astype('float32')

        store = pd.io.pytables.HDFStore(fr)
        self.Fr = store.features.astype('float32')
        store.close()

    def load(self):
        """Load dataset in memory."""
        raise NotImplementedError()


class PredictDataset(Dataset):
    """Test dataset."""

    def __init__(self, fp, fr):
        """
        Constructor.

        Parameters
        ----------
        fp : str
            The name of the HDF5 file containing features for the proteins.

        fr : str
            The name of the HDF5 file containing features for the RNAs.
        """
        super(PredictDataset, self).__init__(fp, fr)

    def load(self):
        """
        Load dataset in memory.

        Return
        ------
        Examples to predict. For each example:
            - p contains the protein features,
            - r contains the RNA features,
            - p_names contains the name of the protein,
            - r_names contains the name of the RNA.

        """
        protein_input_dim = self.Fp.shape[0]
        rna_input_dim = self.Fr.shape[0]
        num_examples = self.Fp.shape[1] * self.Fr.shape[1]
        p = np.zeros((num_examples, protein_input_dim)).astype('float32')
        p_names = []
        r = np.zeros((num_examples, rna_input_dim)).astype('float32')
        r_names = []
        index = 0
        for protein in self.Fp.columns:
            for rna in self.Fr.columns:
                p[index] = self.Fp[protein]
                p_names.append(protein)
                r[index] = self.Fr[rna]
                r_names.append(rna)
                index += 1

        return (p, np.array(p_names), r, np.array(r_names))
