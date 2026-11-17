#!/usr/bin/env python3
"""
Unit-test for the Settings object. The login requests to BioCyc are
mocked, so no data is sent.
"""

import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from cobramod.settings import Settings


class TestBioCycLogin(unittest.TestCase):
    def setUp(self):
        # Settings is a singleton. Keep the credentials of the session, e.g.
        # from credentials.txt, so that other tests can still use them
        self.settings = Settings()
        self.name = self.settings._Settings__biocyc_name
        self.password = self.settings._Settings__biocyc_password
        self.settings._close_biocyc_session()

        self.logins: list[str] = []
        patcher = patch(
            "requests.Session.post",
            side_effect=lambda url, data: self.logins.append(data["email"]),
        )
        patcher.start()
        self.addCleanup(patcher.stop)

    def tearDown(self):
        self.settings._close_biocyc_session()
        self.settings._Settings__biocyc_name = self.name
        self.settings._Settings__biocyc_password = self.password

    def request(self):
        # Every BioCyc request uses this session and logs in if necessary
        return self.settings._biocyc_session

    def test_new_credentials(self):
        # CASE: First request logs in
        self.settings.biocyc_name = "first_user"
        self.settings.biocyc_password = "first_password"
        self.request()
        self.assertEqual(self.logins, ["first_user"])
        self.assertTrue(self.settings.biocyc_logged_in)

        # CASE: Setting the same credentials keeps the session
        self.settings.biocyc_name = "first_user"
        self.settings.biocyc_password = "first_password"
        self.request()
        self.assertEqual(self.logins, ["first_user"])

        # CASE: New credentials close the session and log in again
        self.settings.biocyc_name = "second_user"
        self.assertFalse(self.settings.biocyc_logged_in)
        self.settings.biocyc_password = "second_password"
        self.request()
        self.assertEqual(self.logins, ["first_user", "second_user"])

        # CASE: Only a new password also logs in again
        self.settings.biocyc_password = "third_password"
        self.request()
        self.assertEqual(
            self.logins, ["first_user", "second_user", "second_user"]
        )

    def test_set_biocyc_login(self):
        # CASE: Credentials from a file replace the current ones
        self.settings.biocyc_name = "first_user"
        self.settings.biocyc_password = "first_password"
        self.request()

        with tempfile.TemporaryDirectory() as directory:
            file = Path(directory) / "login.txt"
            file.write_text("file_user\nfile_password\n")
            self.settings.set_biocyc_login(file)

        self.request()
        self.assertEqual(self.logins, ["first_user", "file_user"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
