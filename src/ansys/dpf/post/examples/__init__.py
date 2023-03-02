"""DPF-Post example files."""
import os

# alias files used by the core
from ansys.dpf.core.examples import (
    complex_rst,
    download_all_kinds_of_complexity,
    download_all_kinds_of_complexity_modal,
    download_transient_result,
    electric_therm,
    find_multishells_rst,
    msup_transient,
    multishells_rst,
    simple_bar,
    static_rst,
    steady_therm,
    transient_therm,
)

# must copy files over when using docker
# running_docker = os.environ.get('DPF_DOCKER', False)
