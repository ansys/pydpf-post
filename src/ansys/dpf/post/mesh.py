"""Module containing the ``Mesh`` class."""
from typing import List

from ansys.dpf.core import MeshedRegion

from ansys.dpf.post.data_object import DataObject


class Mesh:
    """Exposes the mesh of the simulation.

    Parameters
    ----------
    meshed_region:
        Instance of :class:`MeshedRegion <ansys.dpf.core.meshed_region.MeshedRegion>` to wrap.

    """

    def __init__(self, meshed_region: MeshedRegion):
        """Initialize this class."""
        self._meshed_region = meshed_region

    @property
    def available_named_selections(self) -> List[str]:
        """Return the available named selection of the mesh."""
        return self._meshed_region.available_named_selections

    @property
    def available_property_fields(self) -> List[str]:
        """Return the available property fields of the mesh."""
        return self._meshed_region.available_property_fields

    @property
    def grid(self):
        """Return the grid of the mesh."""
        return self._meshed_region.grid

    @property
    def nodes(self):
        """Return the nodes of the mesh."""
        return self._meshed_region.nodes

    @property
    def elements(self):
        """Return the elements of the mesh."""
        return self._meshed_region.elements

    @property
    def unit(self):
        """Return the unit of the mesh."""
        return self._meshed_region.unit

    def plot(
        self,
        data: DataObject = None,
        shell_layer: int = None,
        deformation: DataObject = None,
        scale_factor: float = 1.0,
        **kwargs
    ):
        """Plot the mesh.

        Parameters
        ----------
        data :
            Data to plot on the mesh as a contour.
        shell_layer :
            Enum used to set the shell layer if there is data on shell elements.
        deformation :
            Data used to deform the plotted mesh. Must represent a 3D vector field.
        scale_factor :
            Scaling factor applied when deforming the mesh.
        **kwargs : optional
            Additional keyword arguments for the plotter. For additional keyword
            arguments, see ``help(pyvista.plot)``.

        Examples
        --------
        Plot displacement data on a mesh.

        >>> import ansys.dpf.post as dpf
        >>> from ansys.dpf.post import examples
        >>> simulation = dpf.load_simulation(examples.static_rst)
        >>> mesh = simulation.mesh
        >>> displacement = simulation.displacement()
        >>> mesh.plot(displacement)

        Plot the deformed mesh.

        >>> mesh.plot(deformation=displacement, scale_factor=2.0)

        """
        fc = None
        if data:
            fc = data._fc
        deform_by = None
        if deformation:
            deform_by = deformation._fc
        return self._meshed_region.plot(
            field_or_fields_container=fc,
            shell_layers=shell_layer,
            deform_by=deform_by,
            scale_factor=scale_factor,
            **kwargs
        )
