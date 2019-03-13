
import sys
import io
import glob
from setuptools import setup
from distutils.core import Extension
from Cython.Build import cythonize

EXTRA_COMPILE_ARGS = ['-std=c++11']

ext = cythonize([
    Extension('ld_estimator.pairwise',
        extra_compile_args=EXTRA_COMPILE_ARGS,
        sources=['ld_estimator/pairwise.pyx',
            'src/ld.cpp',
            'src/af.cpp',
            'src/confidence.cpp',
            'src/dprime.cpp',
            'src/flip.cpp',
            'src/frequencies.cpp',
            'src/likelihoods.cpp',
            'src/optimize.cpp',
            'src/tallies.cpp',
            'src/utils.cpp',
            ],
        include_dirs=['src/'],
        language='c++')
    ])

setup(name="ld_estimator",
    description='Package for estimating linkage disequilibrium',
    long_description=io.open('README.md', encoding='utf-8').read(),
    long_description_content_type='text/markdown',
    version="1.0.0",
    author="Jeremy McRae",
    author_email="jmcrae@illumina.com",
    url='https://github.com/jeremymcrae/ld_estimator',
    license='MIT',
    packages=['ld_estimator'],
    install_requires=[],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    ext_modules=ext,
    test_suite='tests')
