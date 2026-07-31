|pypi| |python| |downloads| |docs| |bioconda| |galaxy|

CrossMap
=========

CrossMap is a command-line tool for converting genomic coordinates,
annotations, and sequencing alignments between genome assemblies using
chain files.

CrossMap supports many commonly used genomics file formats, including

* BED/BigBed
* GFF/GTF
* VCF/gVCF
* BAM/CRAM
* BigWig
* Wiggle
* BeDGraph
* BAM

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

`https://crossmap.readthedocs.io/ <https://crossmap.readthedocs.io/>`_

Source code
===========

https://github.com/liguowang/CrossMap

Example
=======

Convert a BED file from hg19 to hg38:

::

   CrossMap bed \
      hg19ToHg38.over.chain.gz \
      input.bed \
      output.bed

License
=======

CrossMap is distributed under the GNU General Public License v3.0 or later (GPL-3.0-or-later).

Citation
========

If you use CrossMap in your research, please cite:

Zhao H, Sun Z, Wang J, Huang H, Kocher J-P, Wang L.
CrossMap: a versatile tool for coordinate conversion between genome assemblies.
Bioinformatics. 2014;30(7):1006–1007.

.. |python| image:: https://img.shields.io/pypi/pyversions/CrossMap.svg
   :target: https://pypi.org/project/CrossMap/
   :alt: Python Versions

.. |pypi| image:: https://img.shields.io/pypi/v/CrossMap.svg
   :target: https://pypi.org/project/CrossMap/
   :alt: PyPI Version

.. |bioconda| image:: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg
   :target: https://bioconda.github.io/recipes/crossmap/README.html
   :alt: Install with Bioconda

.. |docs| image:: https://readthedocs.org/projects/crossmap/badge/?version=latest
   :target: https://crossmap.readthedocs.io/
   :alt: Documentation Status

.. |galaxy| image:: https://img.shields.io/badge/Galaxy-use-blue
   :target: https://usegalaxy.eu/root?tool_id=toolshed.g2.bx.psu.edu/repos/iuc/crossmap_bam/crossmap_bam/
   :alt: Available on Galaxy

.. |downloads| image:: https://static.pepy.tech/badge/crossmap
   :target: https://pepy.tech/project/crossmap
   :alt: PyPI Downloads
