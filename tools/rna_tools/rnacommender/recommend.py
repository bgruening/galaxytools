"""Recommend RNAs."""
from __future__ import print_function

import cPickle
import sys
from itertools import izip

from utils import get_serendipity_val

__author__ = "Gianluca Corrado"
__copyright__ = "Copyright 2016, Gianluca Corrado"
__license__ = "MIT"
__maintainer__ = "Gianluca Corrado"
__email__ = "gianluca.corrado@unitn.it"
__status__ = "Production"


class Predictor():
    """Predict interactions."""

    def __init__(self, predict_dataset, trained_model, serendipity_dic=None,
                 output=None):
        """
        Constructor.

        Parameters
        ------
        predict_dataset : data.PredictDataset
            Dataset containing the examples to predict.

        trained_model : str
            File name of the trained model.

        serendipity_dic : dict (default : None)
            Dictionary with serendipy values.

        output : str (default : None)
            Output file. If None then STDOUT.
        """
        self.predict_dataset = predict_dataset
        f = open(trained_model)
        self.model = cPickle.load(f)
        f.close()
        try:
            f = open(serendipity_dic)
            self.serendipity_dic = cPickle.load(f)
            f.close()
        except:
            self.serendipity_dic = None
        self.output = output

    def predict(self):
        """Predict interaction values."""
        # predict the y_hat
        (p, p_names, r, r_names) = self.predict_dataset
        assert p.dtype == 'float32'
        assert r.dtype == 'float32'
        y_hat = self.model.predict(p, r)
        # sort the interactions according to y_hat
        ordering = sorted(range(len(y_hat)),
                          key=lambda x: y_hat[x], reverse=True)
        p_names = p_names[ordering]
        r_names = r_names[ordering]
        y_hat = y_hat[ordering]

        # output to STDOUT
        if self.output is None:
            print("RBP\ttarget\ty_hat\tserendipity")
            if self.serendipity_dic is None:
                for (p_, r_, s_) in izip(p_names, r_names, y_hat):
                    print("%s\t%s\t%.3f\t---" % (p_, r_, s_))
                    sys.stdout.flush()
            else:
                for (p_, r_, s_) in izip(p_names, r_names, y_hat):
                    print("%s\t%s\t%.3f\t%.2f" %
                          (p_, r_, s_,
                           get_serendipity_val(self.serendipity_dic, r_)))
                    sys.stdout.flush()
        # output to file
        else:
            nf = open(self.output, "w")
            nf.write("RBP\ttarget\ty_hat\tserendipity\n")
            if self.serendipity_dic is None:
                for (p_, r_, s_) in izip(p_names, r_names, y_hat):
                    nf.write("%s\t%s\t%.3f\t---\n" % (p_, r_, s_))
            else:
                for (p_, r_, s_) in izip(p_names, r_names, y_hat):
                    nf.write("%s\t%s\t%.3f\t%.2f\n" %
                             (p_, r_, s_,
                              get_serendipity_val(self.serendipity_dic, r_)))
            nf.close()
