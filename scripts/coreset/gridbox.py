"""A docking box around a ligand, and the number of grid points DockThor computes for it.

DockThor reads grid-center, grid-size and grid-rstep as float32 (dockthor/modules/options/options_part.cppm), so
every value here is a float32, written as the shortest decimal that reads back as the same float32. The grid has,
per axis, `static_cast<uint>(size / rstep) + 1` points, the quotient computed in float32
(dockthor/src/Grid.cpp, Grid::Grid), and DockThor stops when the product of the three exceeds --grid-max-points.
"""

import math
import struct
from dataclasses import dataclass
from enum import StrEnum

DEFAULT_PADDING = 8.0  # Å added to each edge of the ligand's bounding box, half on each side
DEFAULT_RSTEP = 0.25  # Å between grid points, DockThor's default --grid-rstep (GRID_DISCRETIZATION_STEP)
DOCKTHOR_MAX_GRID_POINTS = 1_000_000  # DockThor's default --grid-max-points (GRID_NUM_MAX_SIZE)
FLOAT32_MAX_DIGITS = 9  # significant digits that always identify a float32


class BoxShape(StrEnum):
    CUBE = "cube"  # every edge as the longest one
    LIGAND = "ligand"  # each edge as the ligand's own extent on that axis


Vector = tuple[float, float, float]


def float32(value: float) -> float:
    """value rounded to the nearest float32 (as a Python float, which holds it exactly)."""
    return struct.unpack("f", struct.pack("f", value))[0]


def float32_text(value: float) -> float:
    """The float32 of value as the shortest decimal that reads back as that float32: 19.274000549 -> 19.274."""
    target = float32(value)
    for digits in range(1, FLOAT32_MAX_DIGITS + 1):
        shortest = float(f"{target:.{digits}g}")
        if float32(shortest) == target:
            return shortest
    return target


@dataclass(frozen=True)
class GridBox:
    minimum: Vector  # the ligand's bounding box
    maximum: Vector
    center: Vector  # the box written to the TOML, float32
    size: Vector
    rstep: float

    @property
    def extent(self) -> Vector:
        """The edges of the ligand's bounding box."""
        return tuple(high - low for low, high in zip(self.minimum, self.maximum, strict=True))

    @property
    def points(self) -> tuple[int, int, int]:
        """DockThor's grid points per axis."""
        return tuple(grid_points(size, self.rstep) for size in self.size)

    @property
    def total_points(self) -> int:
        return math.prod(self.points)


def grid_points(size: float, rstep: float) -> int:
    """Grid points on an axis of `size` Å, as DockThor computes them: int(size / rstep) + 1 in float32."""
    return int(float32(float32(size) / float32(rstep))) + 1


def grid_box(
    coordinates: list[Vector],
    padding: float = DEFAULT_PADDING,
    shape: BoxShape = BoxShape.CUBE,
    rstep: float = DEFAULT_RSTEP,
) -> GridBox:
    """The box around the atoms: centered on the middle of their bounding box, each edge the bounding box's extent
    plus padding (cube: every edge the longest one), so every atom is at least padding / 2 from the faces."""
    if not coordinates:
        raise ValueError("no atoms")
    minimum = tuple(min(atom[axis] for atom in coordinates) for axis in range(3))
    maximum = tuple(max(atom[axis] for atom in coordinates) for axis in range(3))
    center = tuple(float32_text((low + high) / 2) for low, high in zip(minimum, maximum, strict=True))
    edges = [high - low + padding for low, high in zip(minimum, maximum, strict=True)]
    if shape is BoxShape.CUBE:
        edges = [max(edges)] * 3
    size = tuple(float32_text(edge) for edge in edges)
    return GridBox(minimum, maximum, center, size, float32_text(rstep))
