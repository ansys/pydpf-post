"""Post-processing module.

Using the Data Processing Framework.

Allow to create a result object, then use it to get desired results.

Examples
--------
>>> from ansys.dpf import post
>>> from ansys.dpf.post import examples
>>> solution = post.load_solution(examples.static_rst)
>>> disp = solution.nodal_displacement()

"""

import ansys.dpf.core as core
from ansys.dpf.core.common import locations

from ansys import dpf

try:
    import importlib.metadata as importlib_metadata
except ModuleNotFoundError:  # pragma: no cover
    import importlib_metadata

__version__ = importlib_metadata.version(__name__)


from ansys.dpf.post.common import Grouping as grouping
from ansys.dpf.post.dpf_path import create_path_on_coordinates
from ansys.dpf.post.misc import Report
from ansys.dpf.post.post_utility import load_solution, print_available_keywords

if hasattr(core, "settings") and hasattr(
    core.settings, "set_dynamic_available_results_capability"
):
    core.settings.set_dynamic_available_results_capability(False)
if hasattr(core, "settings") and hasattr(core.settings, "set_default_pyvista_config"):
    core.settings.set_default_pyvista_config()
# dpf.core.start_local_server()
