"""
paratemp: Generic Parallel Tempering (Replica Exchange) Sampler.

Public API:
    - ParallelTempering
    - Replica
    - geometric_temperatures
"""

from .core import (
    ParallelTempering,
    Replica,
    geometric_temperatures,
    BoltzmannDistribution,
    TsallisDistribution,
)

__all__ = [
    "ParallelTempering",
    "Replica",
    "geometric_temperatures",
    "BoltzmannDistribution",
    "TsallisDistribution",
]

__version__ = "0.1.0"
