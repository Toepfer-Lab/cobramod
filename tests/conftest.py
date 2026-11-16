"""
Shared test configuration

The tests download data and add entries to the cross-reference cache in their
data directory. To keep the stored test data unchanged, all test modules work
on a copy of tests/data that is created once per test session.
"""

import shutil
from pathlib import Path

import pytest

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
