---
title: One Click Installation of CTB
layout: post
---

It's done! Thanks to the awesome work of the Galaxy Team, in particular Greg Von Kuster, hours of debugging and 
many features requests we finally have a One Click Installation of ChemicalToolBoX through the Tool Shed.
Search for chemicaltoolbox on the `test toolshed`_ and click install, voila!

To the command line enthusiasts: We have a script for you!

.. code-block:: console

   $ ./scripts/api/install_tools.py --api  <Your Galaxy API Key> \
    -l http://localhost:8080 --url http://testtoolshed.g2.bx.psu.edu/ \
    -o bgruening -r 7e98219aa915  --name chemicaltoolbox \
    --tool-deps --repository-deps \
    --panel-section-name ChemicalToolBoX


.. _`test toolshed`: http://testtoolshed.g2.bx.psu.edu/
