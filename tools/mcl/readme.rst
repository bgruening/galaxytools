====================================
Galaxy wrapper for MCL cluster tools
====================================

The MCL algorithm is short for the Markov Cluster Algorithm, a fast and scalable unsupervised cluster algorithm for graphs (also known as networks) based on simulation of (stochastic) flow in graphs. The algorithm was invented/discovered by Stijn van Dongen (that is, me) at the Centre for Mathematics and Computer Science (also known as CWI) in the Netherlands. The PhD thesis Graph clustering by flow simulation is centered around this algorithm, the main topics being the mathematical theory behind it, its position in cluster analysis and graph clustering, issues concerning scalability, implementation, and benchmarking, and performance criteria for graph clustering in general. The work for this thesis was carried out under supervision of Jan van Eijck and Michiel Hazewinkel. The thesis, technical reports, and preprints can be found in this section. For quickly getting an idea of how MCL operates, consider the flow pictorial at the top of this page, or even better, have a look at an animation of the MCL process. 

This wrapper is copyright 2013 by:
 * Bjoern Gruening
 * Pavan Videm


This prepository contains wrapper for the MCL_ command line tool.

.. _MCL: http://micans.org/mcl/



Stijn van Dongen, Graph Clustering by Flow Simulation, PhD thesis, University of Utrecht, May 2000.
( http://www.library.uu.nl/digiarchief/dip/diss/1895620/inhoud.htm )

Stijn van Dongen, A cluster algorithm for graphs, Technical Report INS-R0010, National Research Institute for Mathematics and Computer Science in the Netherlands, Amsterdam, May 2000.
( http://www.cwi.nl/ftp/CWIreports/INS/INS-R0010.ps.Z )

Stijn van Dongen, Graph clustering via a discrete uncoupling process, Siam Journal on Matrix Analysis and Applications 30-1, p121-141, 2008.

Stijn van Dongen and Cei Abreu-Goodger, Using MCL to extract clusters from networks, in Bacterial Molecular Networks: Methods and Protocols, Methods in Molecular Biology, Vol 804, pages 281—295 (2012). PMID 22144159.



============
Installation
============

MCL will be installed automatically with the Galaxy Tool Shed.



=======
History
=======

interproscan:

 - v0.1: Initial public release of version 12.135


===============================
Wrapper Licence (MIT/BSD style)
===============================

Permission to use, copy, modify, and distribute this software and its
documentation with or without modifications and for any purpose and
without fee is hereby granted, provided that any copyright notices
appear in all copies and that both those copyright notices and this
permission notice appear in supporting documentation, and that the
names of the contributors or copyright holders not be used in
advertising or publicity pertaining to distribution of the software
without specific prior permission.

THE CONTRIBUTORS AND COPYRIGHT HOLDERS OF THIS SOFTWARE DISCLAIM ALL
WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT SHALL THE
CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY SPECIAL, INDIRECT
OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE
OR PERFORMANCE OF THIS SOFTWARE.

