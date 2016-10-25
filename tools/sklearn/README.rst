***************
Galaxy wrapper for scikit-learn library
***************

Contents
========
- `What is scikit-learn?`_
	- `Scikit-learn main package groups`_
	- `Tools offered by this wrapper`_

- `Machine learning workflows`_
	- `Supervised learning workflows`_
	- `Unsupervised learning workflows`_


____________________________


.. _What is scikit-learn?

What is scikit-learn?
===========================

Scikit-learn is an open-source machine learning library for the Python programming language. It offers various algorithms for performing supervised and unsupervised learning as well as data preprocessing and transformation, model selection and evaluation, and dataset utilities. It is built upon SciPy (Scientific Python) library.

Scikit-learn source code can be accessed at https://github.com/scikit-learn/scikit-learn.
Detailed installation instructions can be found at http://scikit-learn.org/stable/install.html


.. _Scikit-learn main package groups:

======
Scikit-learn main package groups
======

Scikit-learn provides the users with several main groups of related operations. 
These are:

- Classification
    - Identifying to which category an object belongs.
- Regression
    - Predicting a continuous-valued attribute associated with an object.
- Clustering
    - Automatic grouping of similar objects into sets.
- Preprocessing
    - Feature extraction and normalization.
- Model selection and evaluation
    - Comparing, validating and choosing parameters and models.
- Dimensionality reduction
    - Reducing the number of random variables to consider.

Each group consists of a number of well-known algorithms from the category. For example, one can find hierarchical, spectral, kmeans, and other clustering methods in sklearn.cluster package.


.. _Tools offered by this wrapper:

===================
Available tools in the current wrapper
===================

The current release of the wrapper offers a subset of the packages from scikit-learn library. You can find:

- A subset of classification metric functions
- Linear and quadratic discriminant classifiers
- Random forest and Ada boost classifiers and regressors
- All the clustering methods
- All supprt vector machine classifiers
- A subset of data preprocessing estimator classes
- Pairwise metric measurement functions

In addition, several tools for performing matrix operations, generating problem-specific datasets, and encoding text and extracting features have been prepared to help the user with more advanced operations.

.. _Machine learning workflows:

Machine learning workflows
===============

Machine learning is about processes. No matter what machine learning algorithm we use, we can apply typical workflows and dataflows to produce more robust models and better predictions.
Here we discuss supervised and unsupervised learning workflows.

.. _Supervised learning workflows:

===================
Supervised machine learning workflows
===================

**What is supervised learning?**

In this machine learning task, given sample data which are labeled, the aim is to build a model which can predict the labels for new observations.
In practice, there are five steps which we can go through to start from raw input data and end up getting reasonable predictions for new samples:

1. Preprocess the data::

    * Change the collected data into the proper format and datatype.
    * Adjust the data quality by filling the missing values, performing 
    required scaling and normalizations, etc.
    * Extract features which are the most meaningfull for the learning task.
    * Split the ready dataset into training and test samples.

2. Choose an algorithm::

    * These factors help one to choose a learning algorithm:
        - Nature of the data (e.g. linear vs. nonlinear data)
        - Structure of the predicted output (e.g. binary vs. multilabel classification)
        - Memory and time usage of the training
        - Predictive accuracy on new data
        - Interpretability of the predictions

3. Choose a validation method
	
	Every machine learning model should be evaluated before being put into practicical use.
	There are numerous performance metrics to evaluate machine learning models. 
	For supervised learning, usually classification or regression metrics are used.

	A validation method helps to evaluate the performance metrics of a trained model in order
	to optimize its performance or ultimately switch to a more efficient model. 
	Cross-validation is a known validation method.

4. Fit a model

   Given the learning algorithm, validation method, and performance metric(s) 
   repeat the following steps:: 

    * Train the model.
    * Evaluate based on metrics. 
    * Optimize unitl satisfied.

5. Use fitted model for prediction::

	This is a final evaluation in which, the optimzed model is used to make predictions
	on unseen (here test) samples. After this, the model is put into production.

.. _Unsupervised learning workflows:

=======================
Unsupervised machine learning workflows
=======================

**What is unsupervised learning?**

Unlike supervised learning and more liklely in real life, here the initial data is not labeled.
The task is to extract the structure from the data and group the samples based on their similarities.
Clustering and dimensionality reduction are two famous examples of unsupervised learning tasks. 

In this case, the workflow is as follows::

    * Preprocess the data (without splitting to train and test).
    * Train a model.
    * Evaluate and tune parameters.
    * Analyse the model and test on real data.

