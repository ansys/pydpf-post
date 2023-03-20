"""This runs at the init of the pytest session

Launch or connect to a persistent local DPF service to be shared in
pytest as a sesson fixture
"""
import os
import re

from ansys.dpf.core.check_version import get_server_version, meets_version
import matplotlib as mpl
import pytest
import pyvista as pv

from ansys.dpf import core
from ansys.dpf.post import examples

# enable off_screen plotting to avoid test interruption
pv.OFF_SCREEN = True
mpl.use("Agg")


def get_lighting():
    """Get lighting configuration.

    Disable lighting when using OSMesa on Windows. See:
    https://github.com/pyvista/pyvista/issues/3185

    """
    pl = pv.Plotter(notebook=False, off_screen=True)
    pl.add_mesh(pv.Sphere())
    pl.show(auto_close=False)
    gpu_info = pl.ren_win.ReportCapabilities()
    pl.close()

    regex = re.compile("OpenGL version string:(.+)\n")
    version = regex.findall(gpu_info)[0]
    return not (os.name == "nt" and "Mesa" in version)


pv.global_theme.lighting = get_lighting()


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
def multishells():
    """Resolve the path of the "multishells.rst" result file."""
    return examples.find_multishells_rst()


@pytest.fixture()
def transient_rst():
    """Resolve the path of the "transient.rst" result file."""
    return examples.download_transient_result()


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


@pytest.fixture()
def modal_frame():
    """Resolve the path of the "rth/rth_electric.rth" result file."""
    return examples.download_modal_frame()


@pytest.fixture(scope="session", autouse=True)
def cleanup(request):
    """Cleanup a testing directory once we are finished."""

    def close_servers():
        core.server.shutdown_all_session_servers()

    request.addfinalizer(close_servers)


SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_6_0 = meets_version(
    get_server_version(core._global_server()), "6.0"
)

SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_5_0 = meets_version(
    get_server_version(core._global_server()), "5.0"
)

# to call at the end
if SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_6_0:
    core.server.shutdown_all_session_servers()
    try:
        core.set_default_server_context(core.AvailableServerContexts.premium)
    except core.errors.DpfVersionNotSupported:
        pass
