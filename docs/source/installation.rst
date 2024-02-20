Requirements
============

CobraMod supports  **Python 3.10 or higher**. Following dependencies are
needed:
 - colorlog
 - cobra>=0.29.0
 - requests>=2.31.0
 - tqdm>=4.62.3
 - openpyxl>=3.1.2
 - webcolors>=1.13
 - pyarrow>=14.0.2

Installation
============

CobraMod can easily be installed using pip. ::

  pip install cobramod

By default, required packages will be installed when installing CobraMod. We
recommend using a virtual environment.

To use the escher visualization tool with cobramod, you can install the optional
dependencies escher-legacy that needs to be used instead of *escher* ::

  pip install cobramod[escher]

.. note::
    For development, install with pip the development branch using the
    argument :code:`-e`. An environment YAML file is provided to use it with a
    conda environment.
