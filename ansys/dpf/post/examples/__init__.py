"""DPF-Post example files"""
import os

# alias files used by the core
from ansys.dpf.core.examples import (
    simple_bar,
    static_rst,
    complex_rst,
    multishells_rst,
    static_rst,
    electric_therm,
    steady_therm,
    transient_therm,
    download_all_kinds_of_complexity_modal,
    download_all_kinds_of_complexity,
    download_transient_result,
    msup_transient,
)

# must copy files over when using docker
# running_docker = os.environ.get('DPF_DOCKER', False)
# breakpoint()
