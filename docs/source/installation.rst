Requirements
============

CobraMod requires  **Python 3.8 or higher** and  the following
packages:

 - cobra >= 0.26.0
 - requests >= 2.24.0
 - webcolors >= 1.11.1
 - Escher >= 1.7.3
 - openpyxl >= 3.0.7
 
Optional:

 - escher-legacy >=1.7.4

Installation
============

CobraMod can easily be installed using pip. ::

  pip install cobramod

By default, required packages will be installed when installing CobraMod. We
recommend using a virtual environment.

To use the escher visualization tool with cobramod, you can install the optional
dependencies. ::

  pip install cobramod[escher]

.. note::
    For development, install with pip the development branch using the
    argument :code:`-e`. An environment YAML file is provided to use it with a
    conda environment.
