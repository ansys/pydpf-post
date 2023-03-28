"""Module containing the ``Mesh`` class.

Mesh
----

"""
from __future__ import annotations

from typing import List

import ansys.dpf.core as dpf

import ansys.dpf.post as post
from ansys.dpf.post import index, locations


class Mesh:
    """Exposes the complete mesh of the simulation."""

    def __init__(self, meshed_region: dpf.MeshedRegion):
        """Initialize this class."""
        self._meshed_region = meshed_region

    def __str__(self):
        """String representation of this class."""
        return str(self._meshed_region).replace("Meshed Region", "Mesh")

    @property
    def available_named_selections(self) -> List[str]:
        """Returns the available named selections of the mesh."""
        return self._meshed_region.available_named_selections

    @property
    def node_ids(self) -> List[int]:
        """Returns the list of node IDs in the mesh."""
        return self._meshed_region.nodes.scoping.ids

    @property
    def element_ids(self) -> List[int]:
        """Returns the list of element IDs in the mesh."""
        return self._meshed_region.elements.scoping.ids

    @property
    def _core_object(self):
        """Returns the underlying PyDPF-Core class:`ansys.dpf.core.MeshedRegion` object."""
        return self._meshed_region

    def plot(self, **kwargs):
        """Plots the Mesh.

        Parameters
        ----------
        kwargs:
            Additional keyword arguments for the plotter. For additional keyword
            arguments, see ``help(pyvista.plot)``.

        Returns
        -------
        A Plotter instance of the current plotting back-end.

        Examples
        --------
        >>> from ansys.dpf import post
        >>> from ansys.dpf.post import examples
        >>> from ansys.dpf.post.common import elemental_properties
        >>> example_path = examples.download_all_kinds_of_complexity()
        >>> simulation = post.StaticMechanicalSimulation(example_path)
        >>> mesh = simulation.mesh
        >>> mesh.plot()

        """
        return self._core_object.plot(**kwargs)

    @property
    def coordinates(self) -> post.DataFrame:
        """Returns the nodal coordinates of the Mesh as a DataFrame.

        Examples
        --------
        >>> from ansys.dpf.post import StaticMechanicalSimulation
        >>> from ansys.dpf.post import examples
        >>> simulation = StaticMechanicalSimulation(examples.static_rst)
        >>> mesh = simulation.mesh
        >>> coord = mesh.coordinates
        >>> print(coord)  # doctest: +NORMALIZE_WHITESPACE
                     results   coord (m)
        node_ids  components
               1           X  1.5000e-02
                           Y  4.5000e-02
                           Z  1.5000e-02
               2           X  1.5000e-02
                           Y  4.5000e-02
                           Z  0.0000e+00
             ...

        """
        from ansys.dpf.post.simulation import vector_component_names

        label = "coord"
        fields_container = dpf.FieldsContainer()
        fields_container.add_field(
            label_space={}, field=self._core_object.nodes.coordinates_field
        )
        return post.DataFrame(
            data=fields_container,
            index=index.MultiIndex(
                indexes=[
                    index.MeshIndex(
                        location=locations.nodal,
                        scoping=self._core_object.nodes.scoping,
                        fc=fields_container,
                    ),
                    index.CompIndex(values=vector_component_names),
                ]
            ),
            columns=index.MultiIndex(
                indexes=[
                    index.ResultsIndex(values=[label], units=fields_container[0].unit)
                ]
            ),
        )
