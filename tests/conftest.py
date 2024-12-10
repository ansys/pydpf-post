"""This runs at the init of the pytest session

Launch or connect to a persistent local DPF service to be shared in
pytest as a sesson fixture
"""
import dataclasses
import os
import pathlib
import re

from ansys.dpf.core.check_version import get_server_version, meets_version
from ansys.dpf.core.examples.downloads import _download_file
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


def save_screenshot(dataframe, suffix=""):
    """Save a screenshot of a dataframe plot, with the current test name."""
    test_path = pathlib.Path(os.environ.get("PYTEST_CURRENT_TEST"))
    dataframe.plot(screenshot=f"{'_'.join(test_path.name.split('::'))}_{suffix}.jpeg")


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
def modalframe():
    """Resolve the path of the "download_modal_frame" result file."""
    return examples.download_modal_frame()


@pytest.fixture()
def simple_cyclic():
    """Resolve the path of the "simple_cyclic.rst" result file."""
    return examples.simple_cyclic


@pytest.fixture()
def multi_stage_cyclic():
    """Resolve the path of the "multi_stage.rst" result file."""
    return examples.download_multi_stage_cyclic_result()


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
def mixed_shell_solid_model():
    """Resolve the path of the "mixed_shell_solid" result file."""
    return _download_file(
        "result_files/extract_shell_layer", "mixed_shell_solid.rst", True, None, False
    )


@pytest.fixture()
def mixed_shell_solid_with_contact_model():
    """Resolve the path of the "mixed_shell_solid_with_contact" result file."""
    return _download_file(
        "result_files/extract_shell_layer",
        "mixed_shell_solid_with_contact.rst",
        True,
        None,
        False,
    )


@pytest.fixture()
def two_cubes_contact_model():
    """Resolve the path of the "two_cubes_contact" result file."""
    return _download_file(
        "result_files/extract_shell_layer", "two_cubes_contact.rst", True, None, False
    )


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
def average_per_body_two_cubes():
    return _download_file(
        "result_files/average_per_body/two_cubes", "file.rst", True, None, False
    )


@pytest.fixture()
def average_per_body_complex_multi_body():
    return _download_file(
        "result_files/average_per_body/complex_multi_body",
        "file.rst",
        True,
        None,
        False,
    )


@dataclasses.dataclass
class ReferenceCsvFilesNodal:
    # reference result with all bodies combined
    # The node ids of nodes at body interfaces are duplicated
    combined: pathlib.Path
    # reference result per body (node ids are unique)
    per_id: dict[str, pathlib.Path]


@dataclasses.dataclass
class ReferenceCsvResult:
    name: str
    has_bodies: bool = True


def get_ref_files(
    root_path: str, n_bodies: int, results: list[ReferenceCsvResult]
) -> dict[str, ReferenceCsvFilesNodal]:
    # Returns a dict of ReferenceCsvFiles for each result_name
    ref_files = {}
    for result in results:
        per_mat_id_dict = {}
        if result.has_bodies:
            for mat in range(1, n_bodies + 1):
                per_mat_id_dict[str(mat)] = _download_file(
                    root_path, f"{result.name}_mat_{mat}.txt", True, None, False
                )
        combined = _download_file(
            root_path, f"{result.name}_combined.txt", True, None, False
        )
        ref_files[result.name] = ReferenceCsvFilesNodal(
            combined=combined, per_id=per_mat_id_dict
        )

    return ref_files


@pytest.fixture()
def average_per_body_complex_multi_body_ref():
    return get_ref_files(
        "result_files/average_per_body/complex_multi_body",
        7,
        results=[ReferenceCsvResult("stress"), ReferenceCsvResult("elastic_strain")],
    )


@pytest.fixture()
def shell_layer_multi_body_ref():
    return get_ref_files(
        "result_files/extract_shell_layer",
        2,
        results=[
            ReferenceCsvResult("stress_top_nodal"),
            ReferenceCsvResult("stress_bot_nodal"),
            ReferenceCsvResult("stress_top_elemental", False),
            ReferenceCsvResult("stress_bot_elemental", False),
        ],
    )


@pytest.fixture()
def average_per_body_two_cubes_ref():
    return get_ref_files(
        "result_files/average_per_body/two_cubes",
        2,
        results=[ReferenceCsvResult("stress"), ReferenceCsvResult("elastic_strain")],
    )


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


@pytest.fixture()
def fluid_fluent_elbow_steady_state():
    """Return paths to fluid fluent mixing elbow steady-state."""
    return examples.download_fluent_mixing_elbow_steady_state()


@pytest.fixture()
def fluid_fluent_elbow_transient():
    """Return paths to fluid fluent mixing elbow transient."""
    return examples.download_fluent_mixing_elbow_transient()


@pytest.fixture()
def fluid_fluent_multi_species():
    """Return paths to fluid fluent multi species."""
    return examples.download_fluent_multi_species()


@pytest.fixture(scope="session", autouse=True)
def cleanup(request):
    """Cleanup a testing directory once we are finished."""

    def close_servers():
        core.server.shutdown_all_session_servers()

    request.addfinalizer(close_servers)


@pytest.fixture
def grpc_server():
    """Starts up a temporary gRPC server"""

    server = core.start_local_server(
        as_global=False, config=core.AvailableServerConfigs.GrpcServer
    )
    yield server
    server.shutdown()


@pytest.fixture(scope="session", autouse=True)
def license_context():
    if SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_6_2:
        with core.LicenseContextManager(
            increment_name="preppost", license_timeout_in_seconds=1.0
        ):
            yield
    else:
        yield


SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_9_1 = meets_version(
    get_server_version(core._global_server()), "9.1"
)

SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_9_0 = meets_version(
    get_server_version(core._global_server()), "9.0"
)

SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_8_0 = meets_version(
    get_server_version(core._global_server()), "8.0"
)

SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_7_1 = meets_version(
    get_server_version(core._global_server()), "7.1"
)

SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_7_0 = meets_version(
    get_server_version(core._global_server()), "7.0"
)

SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_6_2 = meets_version(
    get_server_version(core._global_server()), "6.2"
)

SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_6_0 = meets_version(
    get_server_version(core._global_server()), "6.0"
)

SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_5_0 = meets_version(
    get_server_version(core._global_server()), "5.0"
)

SERVERS_VERSION_GREATER_THAN_OR_EQUAL_TO_4_0 = meets_version(
    get_server_version(core._global_server()), "4.0"
)
