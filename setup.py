from setuptools import setup, find_packages

setup(
    name="cobramod",
    version="0.0.1a0",
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    test_suite="tests",
    python_requires='>=3.7.4',
    install_requires=[
        "cobra>=0.18.1",
        "requests>=2.24.0"],
    url="https://gitlab.com/camborda.s/cobramod"
)
