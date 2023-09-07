"""Module containing the helpers for streamlines.

Streamlines
-----------

"""
from typing import List, Union

from ansys.dpf.core.helpers import streamlines as core_streamlines
from ansys.dpf.core.plotter import DpfPlotter

from ansys.dpf import post


def plot_streamlines(
    dataframe: post.DataFrame,
    sources: List[dict],
    mesh: Union[post.Mesh, None] = None,
    plot_mesh: bool = True,
    mesh_opacity: float = 0.3,
    plot_contour: bool = True,
    contour_opacity: float = 0.3,
    **kwargs,
):
    """Plot streamlines based on a velocity field DataFrame.

    Parameters
    ----------
    dataframe:
        A `post.DataFrame` object containing a velocity field.
    sources:
        A list of dictionaries defining spherical point sources for the streamlines.
        Expected keywords are "center", "radius", and "line_thickness".
    mesh:
        A `post.Mesh` object to use to compute the streamlines instead of the mesh associated with
        the velocity field in the DataFrame.
    plot_contour:
        Whether to plot the field contour along with the streamlines.
    contour_opacity:
        Opacity to use for the field contour in case "plot_contour=True".
    plot_mesh:
        Whether to plot the mesh along the streamlines in case "plot_contour=False".
    mesh_opacity:
        Opacity to use for the mesh in case "plot_contour=False" and "plot_mesh=True".
    **kwargs:

    """
    # Select data to work with
    if mesh:
        meshed_region = mesh._meshed_region
    else:
        meshed_region = dataframe._fc[0].meshed_region
    field = dataframe._fc[0]

    # Initialize the plotter
    plt = DpfPlotter(**kwargs)

    if plot_contour:
        plt.add_field(field=field, opacity=contour_opacity)
    elif plot_mesh:
        plt.add_mesh(meshed_region=meshed_region, opacity=mesh_opacity)

    # Add streamlines for each source
    for source in sources:
        pv_streamline, pv_source = core_streamlines.compute_streamlines(
            meshed_region=meshed_region,
            field=field,
            return_source=True,
            source_radius=source["radius"],
            source_center=source["center"],
        )
        plt.add_streamlines(
            pv_streamline, source=pv_source, radius=source["line_thickness"]
        )

    plt.show_figure(**kwargs)
