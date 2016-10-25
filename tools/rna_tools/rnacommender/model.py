"""Recommender model."""
from __future__ import print_function

import sys

import numpy as np

from theano import function, shared
import theano.tensor as T

__author__ = "Gianluca Corrado"
__copyright__ = "Copyright 2016, Gianluca Corrado"
__license__ = "MIT"
__maintainer__ = "Gianluca Corrado"
__email__ = "gianluca.corrado@unitn.it"
__status__ = "Production"


class Model():
    """Factorization model."""

    def __init__(self, sp, sr, kp, kr, irange=0.01, learning_rate=0.01,
                 lambda_reg=0.01, verbose=True, seed=1234):
        """
        Constructor.

        Parameters
        ----------
        sp : int
            Number of protein features.

        sr : int
            Number of RNA features.

        kp : int
            Size of the protein latent space.

        kr : int
            Size of the RNA latent space.

        irange : float (default : 0.01)
            Initialization range for the model weights.

        learning_rate : float (default : 0.01)
            Learning rate for the weights update.

        lambda_reg : (default : 0.01)
            Lambda parameter for the regularization.

        verbose : bool (default : True)
            Print information at STDOUT.

        seed : int (default : 1234)
            Seed for random number generator.
        """
        if verbose:
            print("Compiling model...", end=' ')
            sys.stdout.flush()

        self.learning_rate = learning_rate
        self.lambda_reg = lambda_reg
        np.random.seed(seed)
        # explictit features for proteins
        fp = T.matrix("Fp", dtype='float32')
        # explictit features for RNAs
        fr = T.matrix("Fr", dtype='float32')
        # Correct label
        y = T.vector("y")

        # projection matrix for proteins
        self.Ap = shared(((.5 - np.random.rand(kp, sp)) *
                          irange).astype('float32'), name="Ap")
        self.bp = shared(((.5 - np.random.rand(kp)) *
                          irange).astype('float32'), name="bp")
        # projection matrix for RNAs
        self.Ar = shared(((.5 - np.random.rand(kr, sr)) *
                          irange).astype('float32'), name="Ar")
        self.br = shared(((.5 - np.random.rand(kr)) *
                          irange).astype('float32'), name="br")
        # generalization matrix
        self.B = shared(((.5 - np.random.rand(kp, kr)) *
                         irange).astype('float32'), name="B")

        # Latent space for proteins
        p = T.nnet.sigmoid(T.dot(fp, self.Ap.T) + self.bp)
        # Latent space for RNAs
        r = T.nnet.sigmoid(T.dot(fr, self.Ar.T) + self.br)
        # Predicted output
        y_hat = T.nnet.sigmoid(T.sum(T.dot(p, self.B) * r, axis=1))

        def _regularization():
            """Normalized Frobenius norm."""
            norm_proteins = self.Ap.norm(2) + self.bp.norm(2)
            norm_rnas = self.Ar.norm(2) + self.br.norm(2)
            norm_b = self.B.norm(2)

            num_proteins = self.Ap.flatten().shape[0] + self.bp.shape[0]
            num_rnas = self.Ar.flatten().shape[0] + self.br.shape[0]
            num_b = self.B.flatten().shape[0]

            return (norm_proteins / num_proteins + norm_rnas / num_rnas +
                    norm_b / num_b) / 3

        # mean squared error
        cost_ = (T.sqr(y - y_hat)).mean()
        reg = lambda_reg * _regularization()
        cost = cost_ + reg

        # compute sgd updates
        g_Ap, g_bp, g_Ar, g_br, g_B = T.grad(
            cost, [self.Ap, self.bp, self.Ar, self.br, self.B])
        updates = ((self.Ap, self.Ap - learning_rate * g_Ap),
                   (self.bp, self.bp - learning_rate * g_bp),
                   (self.Ar, self.Ar - learning_rate * g_Ar),
                   (self.br, self.br - learning_rate * g_br),
                   (self.B, self.B - learning_rate * g_B))

        # training step
        self.train = function(
            inputs=[fp, fr, y],
            outputs=[y_hat, cost],
            updates=updates)
        # test
        self.test = function(
            inputs=[fp, fr, y],
            outputs=[y_hat, cost])

        # predict
        self.predict = function(
            inputs=[fp, fr],
            outputs=y_hat)

        if verbose:
            print("Done.")
            sys.stdout.flush()

    def get_params(self):
        """Return the parameters of the model."""
        return {'Ap': self.Ap.get_value(), 'bp': self.bp.get_value(),
                'Ar': self.Ar.get_value(), 'br': self.br.get_value(),
                'B': self.B.get_value()}
