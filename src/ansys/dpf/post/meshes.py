"""Module containing the ``Meshes`` class.

Meshes
------

"""
from __future__ import annotations

from typing import Union

import ansys.dpf.core as dpf
from ansys.dpf.core import MeshesContainer

from ansys.dpf.post.mesh import Mesh


class Meshes:
    """Container to hold and interact with a split mesh."""

    def __init__(self, meshes_container: MeshesContainer):
        """Initialize this class."""
        self._core_object = meshes_container

    def __getitem__(self, item: Union[int, dict]) -> Mesh:
        """Select a specific mesh based on its position in the container."""
        if isinstance(item, (int, dict)):
            return Mesh(meshed_region=self._core_object.get_mesh(item))
        else:
            ValueError(
                "Access to a specific Mesh of a Meshes requires an index (int) "
                "or a combination of labels (dict)."
            )

    def __str__(self):
        """String representation of this class."""
        return str(self._core_object)

    def select(self, **kwargs) -> Union[Mesh, Meshes, None]:
        """Select a specific mesh based on a combination of property values.

        Parameters
        ----------
        **kwargs:
            A dictionary of property names with associated values to select.
            If values are not defined for a property, all available values are selected.

        Returns
        -------
        A Mesh when the selection results in a unique mesh,
        or a Meshes when several are selected,
        or None when the criteria result in no mesh.

        """
        initial_labels = self._core_object.labels
        # Filter labels
        label_space = {}
        for key in kwargs.keys():
            if key in initial_labels:
                label_space[key] = kwargs[key]
        selected_meshes = self._core_object.get_meshes(label_space=label_space)
        if len(selected_meshes) == 1:
            return Mesh(meshed_region=selected_meshes[0])
        elif len(selected_meshes) == 0:
            return None
        else:
            meshes_container = dpf.MeshesContainer()
            for label in initial_labels:
                meshes_container.add_label(label=label)
            for i, mesh in enumerate(self._core_object):
                if mesh in selected_meshes:
                    new_label_space = self._core_object.get_label_space(i)
                    meshes_container.add_mesh(
                        label_space=new_label_space,
                        mesh=mesh,
                    )
            return Meshes(meshes_container=meshes_container)

    def plot(self, **kwargs):
        """Plots the Meshes."""
        from ansys.dpf.core.plotter import DpfPlotter

        plt = DpfPlotter(**kwargs)
        for mesh in self._core_object:
            plt.add_mesh(meshed_region=mesh, **kwargs)
        return plt.show_figure(**kwargs)
