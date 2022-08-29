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
