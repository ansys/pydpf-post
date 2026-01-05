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

"""Errors specific to DPF-Post."""


class NodalLocationError(ValueError):
    """Provides the error to raise when a location that is set to nodal is invalid."""

    def __init__(self, msg="The location must be nodal."):
        """Initialize this class."""
        ValueError.__init__(self, msg)


class CoreVersionError(ValueError):
    """Provides the error to raise when trying to consume a feature that is unavailable."""

    def __init__(self, version=None, msg=None):
        """Initialize this class."""
        if msg is None:
            msg = """To consume this feature, the ansys-dpf-core
            package version must be later than """
        if version is None:
            version = "0.1.0"
        txt = msg + version
        ValueError.__init__(self, txt)


class LabelSpaceError(ValueError):
    """Provides the error to raise when trying to consume a feature that is unavailable."""

    def __init__(self):
        """Initialize this class."""
        ValueError.__init__(
            self,
            """Arguments display_option and option_id are not correct.
    No corresponding field found to plot.""",
        )


class PandasImportError(ModuleNotFoundError):
    """Error raised when Pandas could not be imported while trying export to a DataFrame."""

    def __init__(
        self,
        msg="To export to a pandas.DataFrame, please install pandas "
        "with :\n pip install pandas",
    ):
        """Initialize this class."""
        ModuleNotFoundError.__init__(self, msg)
