---
title: Installation instruction
layout: default
---

===============
Getting Started
===============

The following instructions are for Linux only.

Prerequisites::

* Python 2.5 or above
* standard C compiler
* Autotools
* CMake


1. Install Mercurial_ at the command prompt if you haven't yet::

	apt-get install mercurial

.. _Mercurial: http://mercurial.selenic.com/

2. Clone the latest `Galaxy platform`_::

	hg clone https://https://bitbucket.org/galaxy/galaxy-central/

.. _Galaxy platform: http://wiki.galaxyproject.org/Admin/Get%20Galaxy

3. Navigate to the galaxy-central folder and update it::
	
	cd ~/galaxy-central
	hg pull
	hg update

4. Create folders for toolshed and dependencies::

	mkdir ~/shed_tools
	mkdir ~/galaxy-central/tool_deps

5. Open universe_wsgi.ini and change the dependencies directory::

	gedit ~/galaxy-central/universe_wsgi.ini


6. Search for ``tool_dependency_dir = None`` and change it to ``tool_dependency_dir = ./tool_deps``

7. (Re-)Start the galaxy daemon::

	GALAXY_RUN_ALL=1 sh run.sh --stop-daemon
	GALAXY_RUN_ALL=1 sh run.sh --daemon

After launching galaxy is accessible via the browser at ``http://localhost:8080/``.

Admin Account
=============

1. Register a new account

2. Promote user to admin
	- open universe_wsgi.ini
	- search ``admin_users = None`` and change it to ``admin_users = YOUR_EMAIL_ADDRESS``

Toolshed
========

1. Installation via webbrowser

- go to the `admin page`_
- select *Search and browse tool sheds*
- Galaxy test tool shed >> Computational chemistry >> chemicaltoolbox
- install chemicaltoolbox

2. Installtion via the shell

- generate an API_ Key
- execute the following code from your terminal

.. code-block:: console

   $ ./scripts/api/install_tools.py --api  <Your Galaxy API Key> \
    -l http://localhost:8080 --url http://testtoolshed.g2.bx.psu.edu/ \
    -o bgruening -r 7e98219aa915  --name chemicaltoolbox \
    --tool-deps --repository-deps \
    --panel-section-name ChemicalToolBoX


.. _admin page: http://localhost:8080/admin
.. _API: http://wiki.galaxyproject.org/Admin/API#Generate_the_Admin_Account_API_Key

VirtualBoX
==========

Dependencies
============


