# Copyright (C) 2020 - 2025 ANSYS, Inc. and/or its affiliates.
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

"""This module contains classes for strain results."""

from ansys.dpf.post.tensor import ComplexTensor, Tensor


class ElasticStrain(Tensor):
    """Defines the elastic strain object, which is a tensor object."""

    def __init__(self, **kwargs):
        """Initialize this class."""
        super().__init__(**kwargs)
        self._operator_name = "EPEL"

    def __str__(self):
        """Return the string representation of this class."""
        return f"{super().__str__()}\nElastic strain object."


class ComplexElasticStrain(ComplexTensor):
    """Defines the complex elastic strain object, which is a tensor object."""

    def __init__(self, **kwargs):
        """Initialize this class."""
        super().__init__(**kwargs)
        self._operator_name = "EPEL"

    def __str__(self):
        """Return the string representation of this class."""
        return f"{super().__str__()}\nComplex elastic strain object."


class PlasticStrain(Tensor):
    """Defines the plastic strain object, which is a tensor object."""

    def __init__(self, **kwargs):
        """Initialize this class."""
        super().__init__(**kwargs)
        self._operator_name = "EPPL"

    def __str__(self):
        """Return the string representation of this class."""
        return f"{super().__str__()}\nPlastic strain object."


class ComplexPlasticStrain(ComplexTensor):
    """Defines the complex plastic strain object, which is a tensor object."""

    def __init__(self, **kwargs):
        """Initialize this class."""
        super().__init__(**kwargs)
        self._operator_name = "EPPL"

    def __str__(self):
        """Return the string representation of this class."""
        return f"{super().__str__()}\nComplex plastic strain object."
