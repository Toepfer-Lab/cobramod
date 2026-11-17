#!/usr/bin/env python3
"""
Unit-test for the Settings object. The login requests to BioCyc are
mocked, so no data is sent.
"""

import tempfile
import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch

from cobramod.settings import Settings


class TestBioCycLogin(unittest.TestCase):
    def setUp(self):
        # Settings is a singleton. Keep the credentials of the session, e.g.
        # from credentials.txt, so that other tests can still use them
        self.settings = Settings()
        self.name = self.settings._Settings__biocyc_name
        self.password = self.settings._Settings__biocyc_password
        self.settings._close_session("BIOCYC")

        # Status code that the mocked BioCyc login returns
        self.status = 200
        self.logins: list[str] = []
        patcher = patch("requests.Session.post", side_effect=self.login)
        patcher.start()
        self.addCleanup(patcher.stop)

    def tearDown(self):
        self.settings._close_session("BIOCYC")
        self.settings._Settings__biocyc_name = self.name
        self.settings._Settings__biocyc_password = self.password

    def login(self, url, data):
        self.logins.append(data["email"])
        response = MagicMock(ok=self.status < 400, status_code=self.status)
        response.json.return_value = {
            "error": "no match for email and password"
        }
        return response

    def request(self):
        # Every BioCyc request uses this session and logs in if necessary
        return self.settings._session("BIOCYC")

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

    def test_rejected_login(self):
        # CASE: Rejected login information is reported and not sent again
        self.status = 401
        self.settings.biocyc_name = "user"
        self.settings.biocyc_password = "wrong_password"

        with self.assertLogs("cobramod.Settings", "WARNING") as logs:
            self.request()
        self.assertIn("no match for email and password", logs.output[0])
        self.assertFalse(self.settings.biocyc_logged_in)
        self.assertEqual(
            self.settings._biocyc_login_error,
            "no match for email and password",
        )

        self.request()
        self.assertEqual(self.logins, ["user"])

        # CASE: New login information is tried again
        self.status = 200
        self.settings.biocyc_password = "correct_password"
        self.assertIsNone(self.settings._biocyc_login_error)
        self.request()
        self.assertEqual(self.logins, ["user", "user"])
        self.assertTrue(self.settings.biocyc_logged_in)

    def test_failed_login(self):
        # CASE: Other errors are tried again with the next request
        self.status = 503
        self.settings.biocyc_name = "user"
        self.settings.biocyc_password = "password"

        with self.assertLogs("cobramod.Settings", "WARNING"):
            self.request()
        self.assertFalse(self.settings.biocyc_logged_in)
        self.assertIsNone(self.settings._biocyc_login_error)

        self.status = 200
        self.request()
        self.assertEqual(self.logins, ["user", "user"])
        self.assertTrue(self.settings.biocyc_logged_in)

    def test_sessions_per_database(self):
        # CASE: Only the BioCyc session logs in
        self.settings.biocyc_name = "user"
        self.settings.biocyc_password = "password"
        kegg = self.settings._session("KEGG")
        self.assertEqual(self.logins, [])

        biocyc = self.request()
        self.assertEqual(self.logins, ["user"])
        self.assertIsNot(kegg, biocyc)
        self.assertIs(self.settings._session("KEGG"), kegg)

        # CASE: New login information only replaces the BioCyc session
        self.settings.biocyc_password = "new_password"
        self.assertIs(self.settings._session("KEGG"), kegg)
        self.assertIsNot(self.request(), biocyc)
        self.settings._close_session("KEGG")


if __name__ == "__main__":
    unittest.main(verbosity=2)
