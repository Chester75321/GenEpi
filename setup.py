# -*- coding: utf-8 -*-
"""
Created on Feb 2018

@author: Chester (Yu-Chuan Chang)
"""

""""""""""""""""""""""""""""""
# import libraries
""""""""""""""""""""""""""""""
from setuptools import setup, find_packages

with open('README.md') as f:
    long_description = f.read()

setup(
    name = 'genepi',
    version = '2.0.8',
    description = 'A package for detecting epsitasis by machine learning',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/Chester75321/GenEpi',
    author = 'Chester (Yu-Chuan Chang)',
    author_email = 'chester75321@gmail.com',
    license = 'MIT',
    project_urls={
        "Source Code": "https://github.com/Chester75321/GenEpi",
        "Documentation": "https://genepi.readthedocs.io/en/latest/",
        "Bug Tracker": "https://github.com/Chester75321/GenEpi/issues",
    },
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Artificial Intelligence',
    ],
    keywords = ['epistasis', 'SNP-SNP interactions', 'GWAS'],
    packages = find_packages(),
    install_requires=[
        'pymysql>=0.8.0',
        'numpy>=1.13.0',
        'scipy>=0.19.0',
        'psutil>=4.3.0',
        'scikit-learn>=0.21.2',
        'matplotlib>=3.1.1',
        'pandas>=0.25.2',
        'seaborn>=0.9.0',
        'PyQt5>=5.14.0',
    ],
    python_requires='>=3',
    include_package_data = True,
    zip_safe = False,
    entry_points = {
        'console_scripts': ['GenEpi=genepi.GenEpi:main'],
    }
)
