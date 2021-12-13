"""This runs at the init of the pytest session

Launch or connect to a persistent local DPF service to be shared in
pytest as a sesson fixture
"""
import os

import pytest
import pyvista as pv

from ansys.dpf import core
from ansys.dpf.post import examples

# enable off_screen plotting to avoid test interruption
pv.OFF_SCREEN = True


# currently running dpf on docker.  Used for testing on CI
running_docker = os.environ.get("DPF_DOCKER", False)


def resolve_test_file(basename, additional_path=""):
    """Resolves a test file's full path based on the base name and the
    environment.

    Normally returns local path unless server is running on docker and
    this repository has been mapped to the docker image at /dpf.
    """
    if running_docker:
        # assumes repository root is mounted at '/dpf'
        test_files_path = "/dpf/tests/testfiles"
        return os.path.join(test_files_path, additional_path, basename)
    else:
        # otherwise, assume file is local
        test_path = os.path.dirname(os.path.abspath(__file__))
        test_files_path = os.path.join(test_path, "testfiles")
        filename = os.path.join(test_files_path, additional_path, basename)
        if not os.path.isfile(filename):
            raise FileNotFoundError(f"Unable to locate {basename} at {test_files_path}")
        return filename


@pytest.fixture()
def allkindofcomplexity():
    """Resolve the path of the "allKindOfComplexity.rst" result file."""
    return examples.download_all_kinds_of_complexity()


@pytest.fixture()
def modalallkindofcomplexity():
    """Resolve the path of the "allKindOfComplexity.rst" result file."""
    return examples.download_all_kinds_of_complexity_modal()


@pytest.fixture()
def simple_bar():
    """Resolve the path of the "ASimpleBar.rst" result file."""
    return examples.simple_bar


@pytest.fixture()
def static_rst():
    """Resolve the path of the "static.rst" result file."""
    return examples.static_rst


@pytest.fixture()
def complex_model():
    """Resolve the path of the "msup/plate1.rst" result file."""
    return examples.complex_rst


@pytest.fixture()
def model_ns():
    """Resolve the path of the "model_with_ns.rst" result file."""
    return examples.multishells_rst


@pytest.fixture()
def plate_msup():
    """Resolve the path of the "msup/plate1.rst" result file.

    Originally:
    UnitTestDataFiles/DataProcessing/expansion/msup/Transient/plate1/file.rst
    """
    return examples.msup_transient


@pytest.fixture()
def rth_transient():
    """Resolve the path of the "rth/rth_transient.rth" result file."""
    return examples.transient_therm


@pytest.fixture()
def rth_steady_state():
    """Resolve the path of the "rth/rth_steady_state.rth" result file."""
    return examples.steady_therm


@pytest.fixture()
def rth_electric():
    """Resolve the path of the "rth/rth_electric.rth" result file."""
    return examples.electric_therm


@pytest.fixture(scope="session", autouse=True)
def cleanup(request):
    """Cleanup a testing directory once we are finished."""

    def close_servers():
        core.server.shutdown_all_session_servers()

    request.addfinalizer(close_servers)
