#!/usr/bin/env python3
"""Setup module

This simple module sets up CobraMod for its installation through pip.
"""
from setuptools import setup, find_packages

setup(
    name="cobramod",
    version="0.2.0",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    package_data={"": ["data/*"]},
    test_suite="tests",
    python_requires=">=3.7.4",
    install_requires=["cobra>=0.18.1", "requests>=2.24.0", "Escher>=1.7.3"],
    url="https://github.com/Toepfer-Lab/cobramod",
)
