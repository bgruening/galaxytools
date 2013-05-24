---
title: New Orphan Packages
layout: post
---

Numpy can be compiled against the atlas_ and lapack_ libraries to gain much better performance.
Today we updated our package for numpy to make use of that.

- http://testtoolshed.g2.bx.psu.edu/view/bgruening/package_atlas_3_10 
- http://testtoolshed.g2.bx.psu.edu/view/bgruening/package_lapack_3_4 
- http://testtoolshed.g2.bx.psu.edu/view/bgruening/package_numpy_1_7 

Additional, the GraphicsMagick_ dependency from OSRA_ was extracted into an own package.

-  http://testtoolshed.g2.bx.psu.edu/view/bgruening/package_graphicsmagick_1_3 


.. _OSRA: http://sourceforge.net/apps/mediawiki/osra
.. _atlas: http://math-atlas.sourceforge.net/
.. _lapack: http://www.netlib.org/lapack/
.. _GraphicsMagick: http://www.graphicsmagick.org/
