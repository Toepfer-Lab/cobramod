#!/usr/bin/env python3
import tempfile
import unittest
from logging import DEBUG
from pathlib import Path
from unittest import TestCase
from unittest.mock import patch

import pandas as pd
from cobra import Metabolite, Reaction
from cobra import __version__ as cobra_version
from numpy import nan

from cobramod import __version__ as cmod_version
from cobramod.core.crossreferences import (
    add2dict_unique,
    add_crossreferences,
    confirms_query,
    get_crossreferences,
    get_reac_prop_with_ec,
    inchikey2pubchem_cid,
    load_cache_from_disk,
    metanetx2ec,
)
from cobramod.debug import debug_log
from cobramod.parsing.db_version import DataVersionConfigurator
from cobramod.settings import Settings

debug_log.setLevel(DEBUG)
data_conf = DataVersionConfigurator()
data_conf.force_same_version = True


class TestCrossReferences(TestCase):
    @patch("requests.get")
    def test_inchikey2pubchem_cid(self, mocked_post):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            mocked_post.return_value.status_code = 200
            mocked_post.return_value.text = "154"

            value = inchikey2pubchem_cid("Test_ID", directory)
            self.assertEqual("154", value)

            # check that the local cache is used
            mocked_post.return_value.text = "155"
            value = inchikey2pubchem_cid("Test_ID", directory)
            self.assertEqual(value, "154")

            value = inchikey2pubchem_cid(["Test_ID2", "Test_ID3"], directory)
            self.assertEqual(["155", "155"], value)

    def test_load_cache_from_disk(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)

            df = load_cache_from_disk("test", directory)
            result = pd.DataFrame({"ID": [], "XRefs": []})

            pd.testing.assert_frame_equal(result, df)

            dict = {"ID": "new_ID", "XRefs": "test"}

            # result = df.append(dict, ignore_index=True)
            result = pd.concat([df, pd.DataFrame([dict])], ignore_index=True)
            (directory / "XRef").mkdir(exist_ok=True)
            result.to_feather(directory / "XRef" / str("test" + ".feather"))
            load_cache_from_disk.cache_clear()
            df = load_cache_from_disk("test", directory)

            pd.testing.assert_frame_equal(result, df)

    @patch("requests.post")
    def test_get_crossreferences(self, mocked_post):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            mocked_post.return_value.status_code = 200
            mocked_post.return_value.json.return_value = {
                "hmdb:HMDB62758": {
                    "InChI": "InChI=1S/C3H8O10P2/"
                    "c4-2(1-12-14(6,7)8)3(5)13-15(9,10)11/"
                    "h2,4H,1H2,(H2,6,7,8)(H2,9,10,11)",
                    "SMILES": "O=C(OP(=O)(O)O)C(O)COP(=O)(O)O",
                    "reference": "chebi:89363",
                    "name": "Glyceric acid 1,3-biphosphate",
                    "InChIkey": "LJQLQCAXBUHEAZ-UHFFFAOYSA-N",
                    "xrefs": [
                        "CHEBI:89363",
                        "chebi:89363",
                        "deprecated:MNXM733151",
                        "hmdb:HMDB0001270",
                        "hmdb:HMDB0062758",
                        "hmdb:HMDB01270",
                        "hmdb:HMDB62758",
                        "sabiork.compound:29",
                        "sabiorkM:29",
                    ],
                    "mnx_id": "MNXM1108074",
                }
            }

            value = get_crossreferences("chem", "hmdb:HMDB62758", directory)

            expected = {
                "inchikey:LJQLQCAXBUHEAZ-UHFFFAOYSA-N",
                "chebi:89363",
                "deprecated:MNXM733151",
                "CHEBI:89363",
                "sabiork.compound:29",
                "sabiorkM:29",
                "hmdb:HMDB62758",
                "hmdb:HMDB0062758",
                "hmdb:HMDB0001270",
                "inchi:InChI=1S/C3H8O10P2/"
                "c4-2(1-12-14(6,7)8)3(5)13-15(9,10)11/"
                "h2,4H,1H2,(H2,6,7,8)(H2,9,10,11)",
                "metanetx.chemical:MNXM1108074",
                "hmdb:HMDB01270",
            }

            self.assertEqual(expected, value)

            mocked_post.return_value.json.return_value = "something_else"
            value = get_crossreferences("chem", "hmdb:HMDB62758", directory)

            self.assertEqual(value, expected)

            value = get_crossreferences(
                "chem", ["hmdb:HMDB62758", "hmdb:HMDB62758"], directory
            )
            self.assertEqual(value, expected)

    def test_confirms_query(self):
        # CASE: Query is the MetaNetX ID or part of a cross-reference
        answer = {
            "mnx_id": "MNXM1364061",
            "xrefs": ["kegg.compound:C00031", "bigg.metabolite:glc__D"],
        }
        self.assertTrue(confirms_query("MNXM1364061", answer))
        self.assertTrue(confirms_query("C00031", answer))
        self.assertTrue(confirms_query("glc__D", answer))
        # CASE: Only exact and case-sensitive matches count
        self.assertFalse(confirms_query("glc__d", answer))
        self.assertFalse(confirms_query("C0003", answer))
        # CASE: Match against the name only (real answer for "m")
        answer = {
            "mnx_id": "MNXM1511782",
            "name": "Petrosaspongiolide m",
            "xrefs": ["hmdb:HMDB0256357"],
        }
        self.assertFalse(confirms_query("m", answer))

    @patch("requests.post")
    def test_get_crossreferences_unprefixed(self, mocked_post):
        # CASE: Answers for queries without prefix are only used if they
        # contain the query as identifier
        mocked_post.return_value.status_code = 200
        mocked_post.return_value.json.return_value = {
            "m": {
                "mnx_id": "MNXM1511782",
                "name": "Petrosaspongiolide m",
                "xrefs": ["hmdb:HMDB0256357"],
            },
            "C00031": {
                "mnx_id": "MNXM1364061",
                "xrefs": ["kegg.compound:C00031", "chebi:4167"],
            },
        }

        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            load_cache_from_disk.cache_clear()

            value = get_crossreferences("chem", ["m", "C00031"], directory)
            self.assertEqual(
                value,
                {
                    "kegg.compound:C00031",
                    "chebi:4167",
                    "metanetx.chemical:MNXM1364061",
                },
            )

            # Rejected answers are cached and not requested again
            mocked_post.reset_mock()
            self.assertEqual(get_crossreferences("chem", "m", directory), set())
            mocked_post.assert_not_called()
            load_cache_from_disk.cache_clear()

    @patch("requests.post")
    def test_get_crossreferences_envipath(self, mocked_post):
        # CASE: enviPath IDs are not sent to MetaNetX
        envipath = (
            "650babc9-9d68-4b73-9332-11972ca26f7b/compound/"
            "98255e64-a8b2-4151-b5f1-08b2ea3a3171"
        )
        mocked_post.return_value.status_code = 200
        mocked_post.return_value.json.return_value = {
            "kegg.compound:C01328": {
                "mnx_id": "MNXM02",
                "xrefs": ["kegg.compound:C01328", f"envipath:{envipath}"],
            }
        }

        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            load_cache_from_disk.cache_clear()

            value = get_crossreferences(
                "chem",
                [
                    "kegg.compound:C01328",
                    f"envipath:{envipath}",
                    f"envipathm:{envipath}",
                ],
                directory,
            )
            self.assertEqual(
                mocked_post.call_args.kwargs["data"]["query_list"],
                "kegg.compound:C01328",
            )
            # enviPath IDs returned as cross-references are kept
            self.assertIn(f"envipath:{envipath}", value)

            mocked_post.reset_mock()
            value = get_crossreferences(
                "chem", f"envipath:{envipath}", directory
            )
            self.assertEqual(value, set())
            mocked_post.assert_not_called()
            load_cache_from_disk.cache_clear()

    @patch("requests.post")
    def test_get_crossreferences_smiles(self, mocked_post):
        # CASE: SMILES are only added if enabled in the settings
        mocked_post.return_value.status_code = 200
        mocked_post.return_value.json.return_value = {
            "hmdb:HMDB62758": {
                "SMILES": "c1:c:c:c:c:c:1",
                "xrefs": ["hmdb:HMDB62758"],
            }
        }
        metabolite = Metabolite()
        metabolite.annotation = {"hmdb": "HMDB62758"}

        settings = Settings()
        settings.add_smiles_as_cross_reference = True
        try:
            with tempfile.TemporaryDirectory() as directory:
                load_cache_from_disk.cache_clear()
                add_crossreferences(metabolite, Path(directory))
        finally:
            settings.add_smiles_as_cross_reference = False
            load_cache_from_disk.cache_clear()

        self.assertEqual(metabolite.annotation["smiles"], "c1:c:c:c:c:c:1")

    @patch("requests.post")
    def test_add_crossreferences_colon_in_identifier(self, mocked_post):
        # CASE: Identifiers can contain ":" (real answer for MNXM102298)
        mocked_post.return_value.status_code = 200
        mocked_post.return_value.json.return_value = {
            "metanetx.chemical:MNXM102298": {
                "mnx_id": "MNXM102298",
                "xrefs": [
                    "metacyc.compound:A-DARA-PEN-1:4",
                    "metacycM:A-DARA-PEN-1:4",
                ],
            }
        }
        metabolite = Metabolite()
        metabolite.annotation = {"metanetx.chemical": "MNXM102298"}

        with tempfile.TemporaryDirectory() as directory:
            load_cache_from_disk.cache_clear()
            add_crossreferences(metabolite, Path(directory))
            load_cache_from_disk.cache_clear()

        self.assertEqual(
            metabolite.annotation["metacyc.compound"], "A-DARA-PEN-1:4"
        )
        self.assertEqual(metabolite.annotation["metacycm"], "A-DARA-PEN-1:4")

    @patch("pandas.read_csv")
    def test_metanetx2ec(self, mock):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            mock.return_value = pd.DataFrame(
                data={
                    "ID": ["test1", "test2"],
                    "mnx_equation": ["test1", "test2"],
                    "reference": ["test1", "test2"],
                    "classifs": ["1.1.1.1", "test2"],
                    "is_balanced": ["test1", "test2"],
                    "is_transport": ["test1", "test2"],
                }
            )

            result = metanetx2ec("test1", directory)
            self.assertEqual("1.1.1.1", result)

            result = metanetx2ec(
                "test1", directory, include_metanetx_specific_ec=True
            )
            self.assertEqual("1.1.1.1", result)

            self.assertRaises(KeyError, metanetx2ec, "test2", directory)

            result = metanetx2ec(
                "test2", directory, include_metanetx_specific_ec=True
            )
            self.assertEqual("test2", result)

    @patch("pandas.read_csv")
    def test_get_reac_prop_with_ec(self, mock):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            mock.return_value = pd.DataFrame(
                data={
                    "ID": ["test1", "test2"],
                    "mnx_equation": ["test1", "test2"],
                    "reference": ["test1", "test2"],
                    "classifs": [nan, "test2"],
                    "is_balanced": ["test1", "test2"],
                    "is_transport": ["test1", "test2"],
                }
            )

            expected = pd.DataFrame(
                data={
                    "ID": ["test2"],
                    "mnx_equation": ["test2"],
                    "reference": ["test2"],
                    "classifs": ["test2"],
                    "is_balanced": ["test2"],
                    "is_transport": ["test2"],
                }
            )

            result = get_reac_prop_with_ec(directory)
            pd.testing.assert_frame_equal(expected, result)

            mock.return_value = "something_else"
            result = get_reac_prop_with_ec(directory)
            pd.testing.assert_frame_equal(expected, result)

    def test_add2dict_unique(self):
        dictonary = {"A": "A", "B": "B"}

        add2dict_unique("A", ["A", "B"], dictonary)
        expected = {"A": ["A", "B"], "B": "B"}

        # check keys
        self.assertCountEqual(expected, dictonary)

        # check values
        for key, value in expected.items():
            self.assertCountEqual(value, dictonary[key])

        add2dict_unique("C", "A", dictonary)
        expected = {"A": ["A", "B"], "B": "B", "C": "A"}

        # check keys
        self.assertCountEqual(expected, dictonary)

        # chack values
        for key, value in expected.items():
            self.assertCountEqual(value, dictonary[key])

    @patch("requests.get")
    @patch("pandas.read_csv")
    @patch("requests.post")
    def test_add_crossreferences(self, mocked_post, mock_pandas, mock_get):
        metabolite = Metabolite()
        metabolite.annotation = {"hmdb": "HMDB62758"}
        reaction = Reaction()
        reaction.annotation = {"hmdb": "HMDB62758"}

        mocked_post.return_value.status_code = 200
        mocked_post.return_value.json.return_value = {
            "hmdb:HMDB62758": {
                "InChI": "InChI=1S/C3H8O10P2/"
                "c4-2(1-12-14(6,7)8)3(5)13-15(9,10)11/"
                "h2,4H,1H2,(H2,6,7,8)(H2,9,10,11)",
                "SMILES": "O=C(OP(=O)(O)O)C(O)COP(=O)(O)O",
                "reference": "chebi:89363",
                "name": "Glyceric acid 1,3-biphosphate",
                "InChIkey": "LJQLQCAXBUHEAZ-UHFFFAOYSA-N",
                "xrefs": [
                    "CHEBI:89363",
                    "chebi:89363",
                    "deprecated:MNXM733151",
                    "hmdb:HMDB0001270",
                    "hmdb:HMDB0062758",
                    "hmdb:HMDB01270",
                    "hmdb:HMDB62758",
                    "sabiork.compound:29",
                    "sabiorkM:29",
                ],
                "mnx_id": "MNXM1108074",
            }
        }

        mock_pandas.return_value = pd.DataFrame(
            data={
                "ID": ["MNXM1108074", "test2"],
                "mnx_equation": ["test1", "test2"],
                "reference": ["test1", "test2"],
                "classifs": ["1.1.1.1", "test2"],
                "is_balanced": ["test1", "test2"],
                "is_transport": ["test1", "test2"],
            }
        )

        mock_get.return_value.text = "154"

        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)

            expected_annotations = {
                "hmdb": [
                    "HMDB01270",
                    "HMDB0001270",
                    "HMDB0062758",
                    "HMDB62758",
                ],
                "metanetx.chemical": "MNXM1108074",
                "sabiorkm": "29",
                "inchi": "InChI=1S/C3H8O10P2/"
                "c4-2(1-12-14(6,7)8)3(5)13-15(9,10)11/"
                "h2,4H,1H2,(H2,6,7,8)(H2,9,10,11)",
                "inchikey": "LJQLQCAXBUHEAZ-UHFFFAOYSA-N",
                "sabiork.compound": "29",
                "chebi": "CHEBI:89363",
                "pubchem.compound": "154",
            }

            add_crossreferences(metabolite, directory)
            # check keys
            self.assertCountEqual(expected_annotations, metabolite.annotation)

            for key, value in expected_annotations.items():
                self.assertCountEqual(value, metabolite.annotation[key])

            add_crossreferences(reaction, directory)
            expected_annotations = {
                "hmdb": [
                    "HMDB0062758",
                    "HMDB62758",
                    "HMDB01270",
                    "HMDB0001270",
                ],
                "metanetx.reaction": "MNXM1108074",
                "sabiorkm": "29",
                "chebi": "CHEBI:89363",
                "metanetx.chemical": "MNXM1108074",
                "sabiork.compound": "29",
                "inchi": "InChI=1S/C3H8O10P2/"
                "c4-2(1-12-14(6,7)8)3(5)13-15(9,10)11/"
                "h2,4H,1H2,(H2,6,7,8)(H2,9,10,11)",
                "inchikey": "LJQLQCAXBUHEAZ-UHFFFAOYSA-N",
                "pubchem.compound": "154",
                "ec-code": "1.1.1.1",
                "brenda": "1.1.1.1",
            }

            # check keys
            self.assertCountEqual(expected_annotations, reaction.annotation)

            for key, value in expected_annotations.items():
                self.assertCountEqual(value, reaction.annotation[key])


if __name__ == "__main__":
    print(f"CobraMod version: {cmod_version}")
    print(f"COBRApy version: {cobra_version}")

    unittest.main(verbosity=2, failfast=True)
