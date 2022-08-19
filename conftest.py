"""This runs at the init of the pytest doctest session

Launch or connect to a persistent local DPF service to be shared in
pytest as a session fixture
"""
import re
import os

import pyvista as pv
import matplotlib as mpl

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
    gpu_info = pl.ren_win.ReportCapabilities();
    pl.close()

    regex = re.compile("OpenGL version string:(.+)\n")
    version = regex.findall(gpu_info)[0]
    return not(os.name == 'nt' and 'Mesa' in version)


pv.global_theme.lighting = get_lighting()
