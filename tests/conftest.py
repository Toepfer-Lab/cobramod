"""
Shared test configuration

The tests download data and add entries to the cross-reference cache in their
data directory. To keep the stored test data unchanged, all test modules work
on a copy of tests/data that is created once per test session.

Login information for the tests can be stored in a ".env" file in the
repository root, see ".env.example".
"""

import os
import shutil
from pathlib import Path

import pytest


def load_env_file(file: Path):
    """
    Sets the variables of a ".env" file as environment variables. Variables
    that are already set, e.g. in the shell, are kept. This must happen
    before CobraMod is imported, as its settings read the login information
    once on import.
    """
    if not file.exists():
        return

    for line in file.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, value = line.split("=", 1)
        value = value.strip().strip("\"'")
        if value:
            os.environ.setdefault(key.strip(), value)


load_env_file(Path(__file__).resolve().parents[1].joinpath(".env"))

dir_source = Path(__file__).resolve().parent.joinpath("data")


@pytest.fixture(scope="session")
def dir_data_copy(tmp_path_factory) -> Path:
    directory = tmp_path_factory.mktemp("cobramod") / "data"
    shutil.copytree(dir_source, directory)
    return directory


@pytest.fixture(scope="module", autouse=True)
def isolated_dir_data(request, dir_data_copy):
    # Test modules define the data directory as module variable 'dir_data'
    if getattr(request.module, "dir_data", None) == dir_source:
        with pytest.MonkeyPatch.context() as monkeypatch:
            monkeypatch.setattr(request.module, "dir_data", dir_data_copy)
            yield
    else:
        yield
