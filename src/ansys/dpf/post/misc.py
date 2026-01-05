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

"""Miscellaneous and report module.

Misc
----

"""
from scooby import Report as ScoobyReport


class Report(ScoobyReport):
    """Generate a report of installed packages for PyDPF-Post."""

    def __init__(self, additional=None, ncol=3, text_width=80, sort=False, gpu=True):
        """Generate a :class:`scooby.Report` instance.

        Parameters
        ----------
        additional : list(ModuleType), list(str), optional
            List of packages or package names to add to the output information.
            The default is ``None``.
        ncol : int, optional
            Number of package-columns in the HTML table. This parameter is only
            used only if ``mode='HTML'`` or ``mode='html'``. The default is ``3``.
        text_width : int, optional
            Text width for non-HTML display modes. The default is ``80``.
        sort : bool, optional
            Whether to sort packages alphabetically. The default is ``False``.
        gpu : bool, optional
            Whether to gather a report on the GPU. The default is ``True``.
            If rendering issues are experienced, set to ``False`` to safely generate
            the report.

        """
        # Mandatory packages.
        core = [
            "pyvista",
            "matplotlib",
            "PIL",
            "pexpect",
            "ansys.dpf.core",
            "ansys.dpf.post",
        ]

        # Optional packages.
        optional = []

        # Information about the GPU - bare except in case there is a rendering
        # bug that the user is trying to report.
        if gpu:
            from pyvista.utilities.errors import GPUInfo

            try:
                extra_meta = [(t[1], t[0]) for t in GPUInfo().get_info()]
            except:
                extra_meta = ("GPU Details", "error")
        else:
            extra_meta = ("GPU Details", "None")

        super().__init__(
            additional=additional,
            core=core,
            optional=optional,
            ncol=ncol,
            text_width=text_width,
            sort=sort,
            extra_meta=extra_meta,
        )


def _connect_any(operator_input, input_value):
    # Workaround to connect any inputs: see
    # https://github.com/ansys/pydpf-core/issues/1670
    operator_input._operator().connect(operator_input._pin, input_value)
