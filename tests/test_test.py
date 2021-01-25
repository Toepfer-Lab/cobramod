"""
This Test checks that the directory data is located in a correct path.
Not all installations will be located at in the directory of cobramod. This
ensure that all installation can use the same data
"""
from pathlib import Path
from unittest import TestCase, main

from cobramod.test import data_dir


class ModuleTest(TestCase):
    def test_directory(self):
        self.assertIn(
            member=str(Path().joinpath("cobramod").joinpath("data")),
            container=str(data_dir),
        )


if __name__ == "__main__":
    main(verbosity=2)
