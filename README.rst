============
ATSAS plugin
============

This plugin provides a wrapper for `ATSAS <http://www.embl-hamburg.de/biosaxs/software.html>`_ - a program suite for small-angle scattering data analysis from biological macromolecules.

.. image:: https://img.shields.io/pypi/v/scipion-em-atsas.svg
        :target: https://pypi.python.org/pypi/scipion-em-atsas
        :alt: PyPI release

.. image:: https://img.shields.io/pypi/l/scipion-em-atsas.svg
        :target: https://pypi.python.org/pypi/scipion-em-atsas
        :alt: License

.. image:: https://img.shields.io/pypi/pyversions/scipion-em-atsas.svg
        :target: https://pypi.python.org/pypi/scipion-em-atsas
        :alt: Supported Python versions

.. image:: https://img.shields.io/sonar/quality_gate/scipion-em_scipion-em-atsas?server=https%3A%2F%2Fsonarcloud.io
        :target: https://sonarcloud.io/dashboard?id=scipion-em_scipion-em-atsas
        :alt: SonarCloud quality gate

.. image:: https://img.shields.io/pypi/dm/scipion-em-atsas
        :target: https://pypi.python.org/pypi/scipion-em-atsas
        :alt: Downloads

Installation
------------

You will need to use 3.0+ version of Scipion to be able to run these protocols. To install the plugin, you have two options:

a) Stable version

.. code-block::

    scipion installp -p scipion-em-atsas

b) Developer's version

    * download repository

    .. code-block::

        git clone -b devel https://github.com/scipion-em/scipion-em-atsas.git

    * install

    .. code-block::

        scipion installp -p /path/to/scipion-em-atsas --devel

**Important:** ATSAS binaries is free for academic users, but you will have to register, download and install them yourself from https://www.embl-hamburg.de/biosaxs/download.html

Configuration variables
-----------------------

*ATSAS_HOME*: Installation path to ATSAS. Default is *software/em/atsas-3.2.0*

Supported versions
------------------

3.0.1, 3.2.0

Protocols
---------

    * convert PDB to SAXS curve

References
----------

    1. Manalastas-Cantos, K., Konarev, P.V., Hajizadeh, N.R., Kikhney, A.G., Petoukhov, M.V., Molodenskiy, D.S., Panjkovich, A., Mertens, H.D.T., Gruzinov, A., Borges, C., Jeffries, C.M., Svergun, D.I., Franke, D. (2021). ATSAS 3.0: expanded functionality and new tools for small-angle scattering data analysis J. Appl. Cryst. 54, 343-355
    2. Svergun, D., Barberato, C., and Koch, M. H. J. (1995). CRYSOL - a Program to Evaluate X-ray Solution Scattering of Biological Macromolecules from Atomic Coordinates. Journal of Applied Crystallography, 28(6): 768–773.
