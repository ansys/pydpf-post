"""Module containing the ``FluidSimulation`` class.

FluidSimulation
---------------

"""
from os import PathLike
from typing import Union

from ansys.dpf import core as dpf
from ansys.dpf.post.simulation import Simulation


class FluidSimulation(Simulation):
    """Base class for fluid type simulations.

    This class provides common methods and properties for all fluid type simulations.
    """

    def __init__(self, result_file: Union[PathLike, str]):
        """Instantiate a mechanical type simulation."""
        model = dpf.Model(result_file)
        data_sources = model.metadata.data_sources
        super().__init__(data_sources=data_sources, model=model)
