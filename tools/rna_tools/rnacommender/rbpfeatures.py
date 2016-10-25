"""Compute the RBP features."""

import re
import sys
import subprocess as sp
import uuid
from os import mkdir
from os import listdir
from os.path import isfile, join
from os import devnull
from shutil import rmtree

import numpy as np

import pandas as pd

import fasta_utils
import pfam_utils

__author__ = "Gianluca Corrado"
__copyright__ = "Copyright 2016, Gianluca Corrado"
__license__ = "MIT"
__maintainer__ = "Gianluca Corrado"
__email__ = "gianluca.corrado@unitn.it"
__status__ = "Production"


class RBPVectorizer():
    """Compute the RBP features."""

    def __init__(self, fasta):
        """
        Constructor.

        Parameters
        ----------
        fasta : str
            Fasta file containing the RBP sequences to predict.
        """
        self.fasta = fasta

        self._mod_fold = "AURA_Human_data/RBP_features/mod"
        self._reference_fisher_scores = \
            "AURA_Human_data/RBP_features/fisher_scores_ref"
        self._train_rbps_file = \
            "AURA_Human_data/RBP_features/rbps_in_train.txt"

        self._temp_fold = "temp_" + str(uuid.uuid4())
        self.pfam_scan = "%s/pfam_scan.txt" % self._temp_fold
        self._dom_fold = "%s/domains" % self._temp_fold
        self._seeds_fold = "%s/seeds" % self._temp_fold
        self._fisher_fold = "%s/fisher_scores" % self._temp_fold

    def _pfam_scan(self):
        """Scan the sequences against the Pfam database."""
        nf = open(self.pfam_scan, "w")
        nf.write(pfam_utils.search_header())

        fasta = fasta_utils.import_fasta(self.fasta)

        if len(fasta) != 1:
            sys.exit("""Fasta file must contain exactly one sequence.""")

        for rbp in sorted(fasta.keys()):
            seq = fasta[rbp]
            text = pfam_utils.sequence_search(rbp, seq)
            nf.write(text)

        nf.close()

    def _overlapping_domains(self):
        """Compute the set of domains contributing to the similarity."""
        reference_domains = set([dom.replace(".mod", "") for dom in
                                 listdir(self._mod_fold) if
                                 isfile(join(self._mod_fold, dom))])

        data = pfam_utils.read_pfam_output(self.pfam_scan)
        if data is None:
            return []

        prot_domains = set([a.split('.')[0] for a in data["hmm_acc"]])
        dom_list = sorted(list(reference_domains & prot_domains))

        return dom_list

    def _prepare_domains(self, dom_list):
        """Select domain subsequences from the entire protein sequences."""
        def prepare_domains(fasta_dic, dom_list, pfam_scan, out_folder):
            out_file_dic = {}
            for acc in dom_list:
                out_file_dic[acc] = open("%s/%s.fa" % (out_folder, acc), "w")

            f = open(pfam_scan)
            f.readline()
            for line in f:
                split = line.split()
                rbp = split[0]
                start = int(split[3])
                stop = int(split[4])
                acc = split[5].split('.')[0]
                if acc in out_file_dic.keys():
                    out_file_dic[acc].write(
                        ">%s:%i-%i\n%s\n" % (rbp, start, stop,
                                             fasta_dic[rbp][start:stop]))
            f.close()

            for acc in dom_list:
                out_file_dic[acc].close()

        mkdir(self._dom_fold)
        fasta = fasta_utils.import_fasta(self.fasta)
        prepare_domains(fasta, dom_list, self.pfam_scan,
                        self._dom_fold)

    def _compute_fisher_scores(self, dom_list):
        """Wrapper for SAM 3.5 get_fisher_scores."""
        def get_fisher_scores(dom_list, mod_fold, dom_fold, fisher_fold):
            for acc in dom_list:
                _FNULL = open(devnull, 'w')
                cmd = "get_fisher_scores run -i %s/%s.mod -db %s/%s.fa" % (
                    mod_fold, acc, dom_fold, acc)
                fisher = sp.check_output(
                    cmd, shell=True, stderr=_FNULL)
                nf = open("%s/%s.txt" % (fisher_fold, acc), "w")
                nf.write(fisher)
                nf.close()

        mkdir(self._fisher_fold)
        get_fisher_scores(dom_list, self._mod_fold, self._dom_fold,
                          self._fisher_fold)

    def _ekm(self, dom_list):
        """Compute the empirical kernel map from the Fisher scores."""
        def process_seg(e):
            """Process segment of a SAM 3.5 get_fisher_scores output file."""
            seg = e.split()
            c = seg[0].split(':')[0]
            m = map(float, seg[3:])
            return c, m

        def read_sam_file(samfile):
            """Read a SAM 3.5 get_fisher_scores output file."""
            f = open(samfile)
            data = f.read()
            f.close()

            columns = []
            m = []
            split = re.split(">A ", data)[1:]
            for e in split:
                c, m_ = process_seg(e)
                columns.append(c)
                m.append(m_)

            m = np.matrix(m)
            df = pd.DataFrame(data=m.T, columns=columns)
            return df

        def dom_features(fisher_fold, dom_list, names=None):
            """Compute the features with respect to a domain type."""
            dfs = []
            for acc in dom_list:
                df = read_sam_file("%s/%s.txt" % (fisher_fold, acc))
                df = df.groupby(df.columns, axis=1).mean()
                dfs.append(df)

            con = pd.concat(dfs, ignore_index=True)

            if names is not None:
                add = sorted(list(set(names) - set(con.columns)))
                fil = sorted(list(set(names) - set(add)))
                con = con[fil]
                for c in add:
                    con[c] = np.zeros(len(con.index), dtype='float64')
                con = con[names]

            con = con.fillna(0.0)
            return con

        f = open(self._train_rbps_file)
        train_rbps = f.read().strip().split('\n')
        f.close()
        ref = dom_features(self._reference_fisher_scores, dom_list,
                           names=train_rbps)
        ekm_ref = ref.T.dot(ref)
        ekm_ref.index = ekm_ref.columns

        sel = dom_features(self._fisher_fold, dom_list)

        ekm_sel = sel.T.dot(sel)
        ekm_sel.index = ekm_sel.columns

        ekm = ref.T.dot(sel)

        for rs in ekm.columns:
            for rr in ekm.index:
                if ekm_ref[rr][rr] > 0 and ekm_sel[rs][rs] > 0:
                    ekm[rs][rr] /= np.sqrt(ekm_ref[rr][rr] * ekm_sel[rs][rs])
        return ekm

    def vectorize(self):
        """Produce the RBP features."""
        # create a temporary folder
        mkdir(self._temp_fold)
        # scan the RBP sequences against Pfam
        self._pfam_scan()
        # determine the accession numbers of the pfam domains needed for
        # computing the features
        dom_list = self._overlapping_domains()
        if len(dom_list) == 0:
            rmtree(self._temp_fold)
            return None
        # prepare fasta file with the sequence of the domains
        self._prepare_domains(dom_list)
        # compute fisher scores using SAM 3.5
        self._compute_fisher_scores(dom_list)
        # compute the empirical kernel map
        ekm = self._ekm(dom_list)
        # remove the temporary folder
        rmtree(self._temp_fold)
        return ekm
