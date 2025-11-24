Settings
*************************

Cobramod's settings object is the central point for controlling its behaviour, setting local save paths, and managing login information.
--------------------
BioCyc-Login
--------------------

To access BioCyc, a license is required. CobraMod offers three ways to provide the necessary login information: a file containing the username and password, environment variables, or passing it directly to the settings object.

File
--------------------
When using a BioCyc login file, it must have two lines. The first line is the username, and the second line is the password. No other text is allowed.

.. code-block:: text

   <Username>
   <Password>

This file can then be used via a settings object, as shown in the following:

.. code-block:: python

   from cobramod import Settings

   settings = Settings()

   settings.SetBioCycLogin("./fileWithLoginInfo.txt")


.. note::
Cobramod also automatically checks the local folder for a file called 'credentials.txt' and, if found, loads the login information from it.
Environment Variables
---------------------

Instead of using a file, CobraMod also supports environment variables. For this, the username needs to be assigned to a variable called 'BIOCYC_NAME', and the password to a variable called 'BIOCYC_PASSWORD'. If no login information is provided, CobraMod checks automatically for the presence of those variables and uses them when necessary.

.. code-block:: text
    BIOCYC_NAME
    BIOCYC_PASSWORD

Via a Settings object
-----------------------

Last but not least, the login information can also be directly passed to the settings object. This can be done via assigning the login information to the ‘biocyc_name’ and ‘biocyc_password’ variables.

.. code-block:: python

   from cobramod import Settings

   settings = Settings()

   settings.biocyc_name = "someone"
   settings.biocyc_password = "APassword"

