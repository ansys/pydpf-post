"""DPF-Post specific errors"""

class NodalLocationError(ValueError):
    """Raised when attempting to """

    def __init__(self, msg="The location must be nodal."):
        ValueError.__init__(self, msg)
