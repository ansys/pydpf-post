"""DPF-Post specific errors"""


class NodalLocationError(ValueError):
    """Raised when attempting to set a location to nodal
    that should not be."""

    def __init__(self, msg="The location must be nodal."):
        ValueError.__init__(self, msg)

class CoreVersionError(ValueError):
    """Raised when attempting to consume a feature that
    is not available in the installed ansys-dpf-core module version."""

    def __init__(self, version=None, msg=None):
        if msg is None:
            msg = """To consume this feature, the ansys-dpf-core
            package version must be higher than """
        if version is None:
            version = '0.1.0'
        txt = msg + version
        ValueError.__init__(self, txt)
