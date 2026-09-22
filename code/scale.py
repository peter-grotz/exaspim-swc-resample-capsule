"""Convert exaSPIM reconstructions between voxel and physical coordinates.

Specimen-space reconstructions are published in voxel units, but resampling has to happen
in microns: the voxel grid is anisotropic, so a fixed step in voxels is a different
physical distance along each axis.

The scale is a property of the acquisition, recorded in ``acquisition.json`` as a
``coordinate_transformations`` entry. The transform stage needs that file for the
registration and carries it forward, so :func:`scale_from_acquisition` reads it directly.

When it is absent, :func:`derive_scale` recovers the same numbers from a matched
voxel/world pair, which this stage always has. That fallback is self-checking: the same
ratio must hold on every axis for every node, and a disagreement means the two
directories are not the same tracing.
"""

from __future__ import annotations

import json
import logging
from collections.abc import Iterator
from pathlib import Path

_LOGGER = logging.getLogger(__name__)

TOLERANCE = 1e-6
"""Relative disagreement tolerated between nodes when deriving the scale."""

SWC_FIELDS = 7
"""Fields in a well-formed SWC node line: id, type, x, y, z, radius, parent."""


class MalformedSwcError(ValueError):
    """Raised when an SWC line does not carry the expected fields."""


class ScaleDerivationError(ValueError):
    """Raised when a voxel-to-physical scale cannot be derived from a matched pair."""


def scale_from_acquisition(acquisition_json: Path) -> tuple[float, float, float] | None:
    """Read the voxel-to-physical scale from an acquisition record.

    Parameters
    ----------
    acquisition_json : Path
        The acquisition carried forward by the transform stage.

    Returns
    -------
    tuple[float, float, float] | None
        Microns per voxel along x, y and z, or ``None`` if the file is missing or holds
        no ``coordinate_transformations`` scale.
    """
    if not acquisition_json.is_file():
        return None
    try:
        payload = json.loads(acquisition_json.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        _LOGGER.warning("Could not read %s: %s", acquisition_json, error)
        return None

    for transform in _iter_transforms(payload):
        if transform.get("type") == "scale" and len(transform.get("scale", ())) == 3:
            scale = tuple(float(value) for value in transform["scale"])
            _LOGGER.info("Read voxel scale %s um from %s", scale, acquisition_json.name)
            return scale
    _LOGGER.warning("No coordinate_transformations scale in %s", acquisition_json.name)
    return None


def _iter_transforms(payload: object) -> "Iterator[dict]":
    """Yield every ``coordinate_transformations`` entry in a nested payload.

    The entries live per-tile inside the acquisition, so the structure is walked rather
    than indexed at a fixed depth.

    Parameters
    ----------
    payload : object
        Decoded acquisition JSON, at any depth.

    Yields
    ------
    dict
        Each transform entry found.
    """
    if isinstance(payload, dict):
        for key, value in payload.items():
            if key == "coordinate_transformations" and isinstance(value, list):
                yield from (item for item in value if isinstance(item, dict))
            else:
                yield from _iter_transforms(value)
    elif isinstance(payload, list):
        for item in payload:
            yield from _iter_transforms(item)


def _node_fields(line: str, path: Path) -> list[str]:
    """Split an SWC node line, rejecting one that is truncated.

    A partially written file would otherwise fail with an opaque ``IndexError`` deep in
    the coordinate handling.

    Parameters
    ----------
    line : str
        The line to split.
    path : Path
        File the line came from, for the error message.

    Returns
    -------
    list[str]
        The line's whitespace-separated fields.

    Raises
    ------
    MalformedSwcError
        If the line has fewer than :data:`SWC_FIELDS` fields.
    """
    fields = line.split()
    if len(fields) < SWC_FIELDS:
        raise MalformedSwcError(
            f"{path.name} has a node line with {len(fields)} fields, expected "
            f"{SWC_FIELDS}: {line.strip()[:60]!r}"
        )
    return fields


def _read_points(path: Path, limit: int) -> list[tuple[float, float, float]]:
    """Read the leading coordinates of an SWC.

    Parameters
    ----------
    path : Path
        SWC to read.
    limit : int
        Maximum number of nodes to read.

    Returns
    -------
    list[tuple[float, float, float]]
        The ``(x, y, z)`` of each node read, in file order.
    """
    points: list[tuple[float, float, float]] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = _node_fields(line, path)
            points.append((float(fields[2]), float(fields[3]), float(fields[4])))
            if len(points) >= limit:
                break
    return points


def derive_scale(
    voxel_swc: Path,
    world_swc: Path,
    sample_nodes: int = 200,
) -> tuple[float, float, float]:
    """Derive the voxel-to-physical scale from a matched pair of reconstructions.

    Parameters
    ----------
    voxel_swc : Path
        The reconstruction in voxel coordinates.
    world_swc : Path
        The same reconstruction in physical coordinates.
    sample_nodes : int, optional
        How many leading nodes to compare, by default 200.

    Returns
    -------
    tuple[float, float, float]
        Microns per voxel along x, y and z.

    Raises
    ------
    ScaleDerivationError
        If either file is empty, the two disagree in length, or the ratio is not constant
        across the sampled nodes — any of which means they are not the same tracing.
    """
    voxel = _read_points(voxel_swc, sample_nodes)
    world = _read_points(world_swc, sample_nodes)
    if not voxel or not world:
        raise ScaleDerivationError(f"No nodes read from {voxel_swc.name} or {world_swc.name}")
    if len(voxel) != len(world):
        raise ScaleDerivationError(
            f"{voxel_swc.name} and {world_swc.name} differ in length "
            f"({len(voxel)} vs {len(world)}); they are not the same tracing"
        )

    scale: list[float] = []
    for axis in range(3):
        ratios = [w[axis] / v[axis] for v, w in zip(voxel, world) if abs(v[axis]) > TOLERANCE]
        if not ratios:
            raise ScaleDerivationError(
                f"Axis {axis} is zero throughout {voxel_swc.name}; cannot derive a scale"
            )
        first = ratios[0]
        if any(abs(ratio - first) > TOLERANCE * max(1.0, abs(first)) for ratio in ratios):
            raise ScaleDerivationError(
                f"Axis {axis} scale is not constant between {voxel_swc.name} and "
                f"{world_swc.name}; they are not the same tracing"
            )
        scale.append(first)

    result = (scale[0], scale[1], scale[2])
    _LOGGER.info("Derived voxel scale %s um from %s", result, voxel_swc.name)
    return result


def scale_swc(source: Path, destination: Path, factors: tuple[float, float, float]) -> int:
    """Write ``source`` with each coordinate multiplied by ``factors``.

    Comment lines are preserved, so an ``# OFFSET`` header survives unchanged. Note that
    such a header is in the *source* units and is not rescaled.

    Parameters
    ----------
    source : Path
        SWC to read.
    destination : Path
        SWC to write. Parent directories are created.
    factors : tuple[float, float, float]
        Multipliers for x, y and z.

    Returns
    -------
    int
        Number of nodes written.
    """
    destination.parent.mkdir(parents=True, exist_ok=True)
    written = 0
    with source.open("r", encoding="utf-8") as reader, destination.open(
        "w", encoding="utf-8"
    ) as writer:
        for line in reader:
            if line.startswith("#") or not line.strip():
                writer.write(line)
                continue
            fields = _node_fields(line, source)
            for axis in range(3):
                fields[2 + axis] = f"{float(fields[2 + axis]) * factors[axis]:.6f}"
            writer.write(" ".join(fields) + "\n")
            written += 1
    return written
