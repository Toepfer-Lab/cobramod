#!/usr/bin/env python3
from cobramod.debug import debug_log
from pathlib import Path
import json
import requests


def _get_json_bigg(
    directory: Path, identifier: str, model_id: str, object_type: str
) -> dict:
    """
    Searchs in given parent directory if data is located in their respective
    database directory. If not, data will be retrievied from the corresponding
    database. Returns root of given identifier.

    Args:
        directory (Path): Path to directory where data is located.
        identifier (str): identifier for given database.
        model_id (str): Name of model to get the information
        object_type (str): Either 'reaction' or 'metabolite'.

    Raes:
        Warning: If object is not available in given database
        NotADirectoryError: If parent directory is not found

    Returns:
        ET.Element: root of XML file
    """
    if directory.exists():
        data_dir = directory.joinpath("BIGG")
        if not data_dir.exists():
            data_dir.mkdir()
        filename = data_dir.joinpath(f"{identifier}.json")
        debug_log.debug(
            f'Searching {object_type} "{identifier}" in directory "BIGG"'
        )
        try:
            with open(file=filename, mode="r") as f:
                unformatted_data = f.read()
            debug_log.debug(f"Identifier '{identifier}' found.")
            return json.loads(s=unformatted_data)
        except FileNotFoundError:
            debug_log.debug(
                f'{object_type.capitalize()} "{identifier}" not found in '
                f'directory "BIGG".'
            )
            # Retrieve from URL
            url_text = (
                f"http://bigg.ucsd.edu/api/v2/models/{model_id}/{object_type}s"
                f"/{identifier}"
            )
            debug_log.debug(f"Searching in {url_text}")
            r = requests.get(url_text)
            if r.status_code == 404:
                msg = f'"{identifier}" not found as {object_type}'
                debug_log.error(msg)
                raise Warning(msg)
            else:
                unformatted_data = r.text
                debug_log.info(
                    f'Object "{identifier}" found. Saving in '
                    f'directory "BIGG".'
                )
                with open(file=filename, mode="w+") as f:
                    f.write(unformatted_data)
                return json.loads(s=unformatted_data)
    else:
        msg = "Directory not found"
        debug_log.critical(msg)
        raise NotADirectoryError(msg)
