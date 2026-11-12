"""coreset.gridbox: the box around a ligand and DockThor's grid-point count."""

import pytest

from coreset.gridbox import BoxShape, float32, float32_text, grid_box, grid_points

# A ligand whose bounding box is 10 x 20 x 30 Å, from (0, 0, 0) to (10, 20, 30).
ATOMS = [(0.0, 0.0, 0.0), (10.0, 5.0, 1.0), (3.0, 20.0, 30.0)]


def test_cube_takes_the_longest_edge_plus_padding():
    box = grid_box(ATOMS, padding=8.0)
    assert box.extent == (10.0, 20.0, 30.0)
    assert box.size == (38.0, 38.0, 38.0)
    assert box.center == (5.0, 10.0, 15.0)  # the middle of the bounding box, not the mean of the atoms


def test_ligand_shape_keeps_each_edge():
    box = grid_box(ATOMS, padding=8.0, shape=BoxShape.LIGAND)
    assert box.size == (18.0, 28.0, 38.0)


def test_every_atom_keeps_half_the_padding_to_the_faces():
    atoms = [(1.5, -2.25, 7.0), (4.125, 3.0, 9.5), (2.0, 0.5, 8.0), (1.6, 2.9, 7.2)]
    for shape in BoxShape:
        box = grid_box(atoms, padding=6.0, shape=shape)
        for atom in atoms:
            for axis in range(3):
                margin = box.size[axis] / 2 - abs(atom[axis] - box.center[axis])
                assert margin >= 3.0 - 1e-5


@pytest.mark.parametrize(
    ("size", "rstep", "points"),
    [
        (38.0, 0.25, 153),  # an exact multiple: 152 intervals
        (19.274, 0.25, 78),  # 77.096 intervals, truncated
        (28.585, 0.25, 115),  # 3ag9 with 8 Å of padding
        (26.585, 0.2987, 90),  # 3ag9 in the current set: 89 intervals
        (0.1, 0.25, 1),
    ],
)
def test_grid_points_per_axis(size, rstep, points):
    assert grid_points(size, rstep) == points


def test_total_points_over_dockthors_limit():
    box = grid_box([(0.0, 0.0, 0.0), (20.585, 1.0, 1.0)], padding=8.0)
    assert box.points == (115, 115, 115)
    assert box.total_points == 1_520_875


def test_float32_quotient_not_float64():
    # 0.3 / 0.1 is 2.9999999999999996 in float64 (2 intervals) but exactly 3.0 in float32, as DockThor computes it.
    assert int(0.3 / 0.1) == 2
    assert grid_points(0.3, 0.1) == 4


def test_float32_text_is_the_shortest_round_trip():
    assert float32_text(19.274000549316406) == 19.274
    assert float32_text(26.5055) == 26.5055
    assert float32_text(0.1) == 0.1
    for value in (13.7395, 1e-3, 123456.78, -46.043, 1 / 3):
        text = float32_text(value)
        assert float32(text) == float32(value)
        assert len(repr(text)) <= len(repr(float32(value)))


def test_no_atoms():
    with pytest.raises(ValueError):
        grid_box([])
