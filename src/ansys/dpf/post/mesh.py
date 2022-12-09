"""Module containing the ``Mesh`` class.
"""
from ansys.dpf.core import MeshedRegion


class Mesh:
    def __init__(self, meshed_region: MeshedRegion):
        self._meshed_region = meshed_region
