#!/usr/bin/env python3
"""Setup module

This simple module sets up CobraMod for its installation through pip.
"""
from setuptools import setup, find_packages

with open(file="README.md", mode="r") as f:
    readme = "".join(f.readlines())

setup(
    name="cobramod",
    version="1.0.2",
    description="Python package for pathway-centric modification and extension"
    + " of genome-scale metabolic networks",
    long_description=readme,
    long_description_content_type="text/markdown",
    author="Stefano Camborda La Cruz, "
    + "Jan-Niklas Weder, "
    + "Nadine Töpfer",
    author_email="toepfer@ipk-gatersleben.de",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    package_data={"": ["data/*"]},
    test_suite="tests",
    python_requires=">=3.7.4",
    install_requires=[
        "cobra>=0.18.1",
        "requests>=2.24.0",
        "Escher>=1.7.3",
        "openpyxl>=3.0.7",
        "webcolors>=1.11.1",
        "tqdm>=4.62.3",
        "pyarrow>=6.0.1",
    ],
    url="https://github.com/Toepfer-Lab/cobramod",
    project_urls={
        "Documentation": "https://cobramod.readthedocs.io/",
        "Bug Tracker": "https://github.com/Toepfer-Lab/cobramod/issues",
    },
    keywords=[
        "genome-scale metabolic model",
        "constraint-based modelling",
        "COBRApy",
        "Escher",
        "metabolic model curation",
    ],
    classifiers=[
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    license="GPL v.3.0",
    platforms=[""],
)
