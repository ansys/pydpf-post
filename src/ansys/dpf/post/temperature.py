# Copyright (C) 2020 - 2026 ANSYS, Inc. and/or its affiliates.
# SPDX-License-Identifier: MIT
#
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

"""This module contains classes for temperature results ."""

from ansys.dpf.post.common import _AvailableKeywords
from ansys.dpf.post.scalar import ComplexScalar, Scalar
from ansys.dpf.post.vector import Vector


class StructuralTemperature(Scalar):
    """Defines the structural temperature object, which is a scalar object."""

    def __init__(self, **kwargs):
        """Initialize this class."""
        super().__init__(**kwargs)
        self._operator_name = "BFE"

    def __str__(self):
        """Return the string representation of this class."""
        return f"{super().__str__()}\nStructural temperature object."


class ComplexStructuralTemperature(ComplexScalar):
    """Defines the complex structural temperature object, which is a scalar object."""

    def __init__(self, **kwargs):
        """Initialize this class."""
        super().__init__(**kwargs)
        self._operator_name = "BFE"

    def __str__(self):
        """Return the string representation of this class."""
        return f"{super().__str__()}\nComplex temperature object."


class Temperature(Scalar):
    """Defines the temperature object, which is a scalar object, for thermal analysis."""

    def __init__(self, **kwargs):
        """Initialize this class."""
        super().__init__(**kwargs)
        self._operator_name = "TEMP"

        # disable element scoping
        if _AvailableKeywords.element_scoping in kwargs:
            raise Exception(
                "Element scoping is not available with thermal/electric results."
            )
        self.definition._Definition__element_scoping_locked = True

    def __str__(self):
        """Return the string representation of this class."""
        return f"{super().__str__()}\nTemperature object."


class HeatFlux(Vector):
    """Defines the heat flux object, which is a scalar object, for thermal analysis."""

    def __init__(self, **kwargs):
        """Initialize this class."""
        super().__init__(**kwargs)
        self._operator_name = "TF"

        # disable element scoping
        if _AvailableKeywords.element_scoping in kwargs:
            raise Exception(
                "Element scoping is not available with thermal/electric results."
            )
        self.definition._Definition__element_scoping_locked = True

    def __str__(self):
        """Return the string representation of this class."""
        return f"{super().__str__()}\nHeat flux object."
