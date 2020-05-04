import sys, os, platform, glob
from distutils.core import setup
from setuptools import *

"""
Setup script for CrossMap  -- Lift over genomics coordinates between assemblies
"""


def main():
    setup(  name = "CrossMap",
            version = "0.4.2",
            python_requires='>=3.5',
            py_modules = [ 'psyco_full' ],
            packages = find_packages( 'lib' ),
            package_dir = { '': 'lib' },
            package_data = { '': ['*.ps'] },
            scripts = glob.glob( "bin/*.py"),
            ext_modules = [],
            test_suite = 'nose.collector',
            setup_requires = ['nose>=0.10.4','cython>=0.12'],
            author = "Liguo Wang",
            author_email ="wangliguo78@gmail.com",
            platforms = ['Linux','MacOS'],
            requires = ['cython (>=0.17)'],
            install_requires = ['cython>=0.17','pysam','bx-python','pyBigWig'], 
            description = " Lift over genomics coordinates between assemblies",
            url = "http://crossmap.sourceforge.net/",
            zip_safe = False,
            dependency_links = [],
            classifiers=[
                'Development Status :: 5 - Production/Stable',
                'Environment :: Console',
                'Intended Audience :: Science/Research',
                'License :: OSI Approved :: GNU General Public License v2 (GPLv2)',
                'Operating System :: MacOS :: MacOS X',
                'Operating System :: POSIX',
                'Programming Language :: Python',
                'Topic :: Scientific/Engineering :: Bio-Informatics',
            ],
            
            keywords='Genome coordinates lift over',
             )


if __name__ == "__main__":
    main()
