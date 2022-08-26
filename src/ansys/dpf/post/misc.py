"""Miscellaneous and report module."""
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
