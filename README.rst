|docs| |bioconda| |galaxy|

CrossMap
=========

CrossMap is a command-line tool for converting genomic coordinates,
annotation files, and sequencing alignments between genome assemblies
using UCSC chain files.

CrossMap supports many commonly used genomics file formats, including

* BED
* GFF/GTF
* VCF
* BAM/CRAM
* BigWig
* Wiggle
* BEDGraph

Installation
============

Install the latest release from PyPI:

::

    pip install CrossMap

Or install the latest development version from GitHub:

::

    pip install git+https://github.com/liguowang/CrossMap.git

Documentation
=============

https://crossmap.readthedocs.io/

Source code
===========

https://github.com/liguowang/CrossMap

Example
=======

Convert a BED file from hg19 to hg38:

::

    CrossMap bed hg19ToHg38.over.chain.gz input.bed output.bed

License
=======

CrossMap is distributed under the GNU General Public License v3.0 or later (GPL-3.0-or-later).

.. |docs| image:: https://readthedocs.org/projects/crossmap/badge/?version=latest
   :target: https://crossmap.readthedocs.io/
   :alt: Documentation Status

.. |bioconda| image:: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg
   :target: https://bioconda.github.io/recipes/crossmap/README.html
   :alt: Install with Bioconda

.. |galaxy| image:: https://img.shields.io/badge/Galaxy-use-blue
   :target: https://usegalaxy.eu/root?tool_id=toolshed.g2.bx.psu.edu/repos/iuc/crossmap_bam/crossmap_bam/
   :alt: Available on Galaxy
