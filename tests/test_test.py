from pathlib import Path
import unittest

from cobramod.test import data_dir


class TestTest(unittest.TestCase):
    def test_directory(self):
        self.assertEqual(
            first=str(data_dir),
            second=str(
                Path.cwd()
                .joinpath("src")
                .joinpath("cobramod")
                .joinpath("data")
            ),
        )


if __name__ == "__main__":
    unittest.main(verbosity=2)
