
import io
from setuptools import setup

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
    test_suite='tests')
