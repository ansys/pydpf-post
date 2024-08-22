from enum import Enum


class ResultCategory(Enum):
    """Enum for available result categories."""

    scalar = 1
    vector = 2
    matrix = 3
    principal = 4
    equivalent = 5
