"""Resample exaSPIM reconstructions and annotate the CCF-space ones.

Two resampling passes, both at 10 um:

* **CCF space** — the aligned reconstructions, then annotated with CCF structure ids and
  written as MouseLight JSON.
* **Specimen space** — resampled in *physical* coordinates and converted back to voxels.
  Specimen coordinates are voxels on an anisotropic grid, so resampling them directly
  would space nodes differently along each axis; see :mod:`scale`.

Code Ocean glue around ``neuron-tracing-utils`` and ``aind-morphology-utils``. Their
pinned numeric stack keeps them out of the shared library, which supplies the stage
metadata record.
"""

import argparse
import json
import logging
import os
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

from aind_morphology_utils.ccf_annotation import CCFMorphologyMapper
from aind_morphology_utils.utils import read_swc
from aind_morphology_utils.writers import MouseLightJsonWriter
from exaspim_swc_processing.naming import ReconstructionNameError, parse_stem
from exaspim_swc_processing.stage import build_stage_process, resolve_code, write_stage_process
from scale import derive_scale, scale_from_acquisition, scale_swc

DATA_DIR = Path(os.environ.get("DATA_DIR", "/data"))
RESULTS_DIR = Path(os.environ.get("RESULTS_DIR", "/results"))
SCRATCH_DIR = Path(os.environ.get("SCRATCH_DIR", "/scratch"))
STEP_NAME = "aligned_swc_processing"
DEFAULT_SPACING_UM = 10.0

logger = logging.getLogger(STEP_NAME)


def parse_args() -> argparse.Namespace:
    """Read the App Builder parameters.

    Returns
    -------
    argparse.Namespace
        Parsed arguments.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--spacing-um",
        type=float,
        default=float(os.environ.get("NODE_SPACING_UM", DEFAULT_SPACING_UM)),
        help="Node spacing in microns, applied to both coordinate spaces.",
    )
    parser.add_argument(
        "--ccf-resolution-um",
        type=float,
        default=float(os.environ.get("CCF_RESOLUTION_UM", DEFAULT_SPACING_UM)),
    )
    parser.add_argument("--experimenters", default=os.environ.get("EXPERIMENTERS", ""))
    parser.add_argument("--fail-fast", action="store_true")
    return parser.parse_args()


def resample(source: Path, destination: Path, spacing_um: float) -> None:
    """Resample every reconstruction in a directory to a fixed node spacing.

    Parameters
    ----------
    source : Path
        Directory of reconstructions, in microns.
    destination : Path
        Directory to write to.
    spacing_um : float
        Target spacing between nodes.
    """
    destination.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [
            sys.executable,
            "-m",
            "neuron_tracing_utils.resample",
            "--input",
            str(source),
            "--output",
            str(destination),
            "--spacing",
            str(spacing_um),
        ],
        check=True,
    )


def find_dir(*candidates: str) -> Path | None:
    """Return the first candidate that exists under the mounted inputs.

    Parameters
    ----------
    *candidates : str
        Paths relative to the data directory, in priority order.

    Returns
    -------
    Path | None
        The first that is a directory, or ``None``.
    """
    for candidate in candidates:
        path = DATA_DIR / candidate
        if path.is_dir():
            return path
    return None


def matched_pair(voxel_dir: Path, world_dir: Path) -> tuple[Path, Path]:
    """Find one reconstruction present in both coordinate spaces.

    Parameters
    ----------
    voxel_dir : Path
        Directory of voxel-space reconstructions.
    world_dir : Path
        Directory of physical-space reconstructions.

    Returns
    -------
    tuple[Path, Path]
        The voxel and physical forms of the same reconstruction.

    Raises
    ------
    FileNotFoundError
        If no stem appears in both directories.
    """
    for voxel in sorted(voxel_dir.glob("*.swc")):
        world = world_dir / voxel.name
        if world.is_file():
            return voxel, world
    raise FileNotFoundError(
        f"No reconstruction appears in both {voxel_dir} and {world_dir}; "
        "cannot derive the voxel scale"
    )


def resample_specimen_space(spacing_um: float, output_dir: Path) -> dict[str, object]:
    """Resample specimen-space reconstructions at a physical spacing.

    Resampling runs on the physical-coordinate forms so the spacing is isotropic in
    microns, then the result is converted back to voxels for publication.

    Parameters
    ----------
    spacing_um : float
        Target spacing between nodes, in microns.
    output_dir : Path
        Directory to write voxel-space resampled reconstructions to.

    Returns
    -------
    dict[str, object]
        Facts for the stage record: the scale used and how many were written.
    """
    voxel_dir = find_dir("refinement/final-voxel", "swc_refinement/final-voxel", "final-voxel")
    world_dir = find_dir("refinement/final-world", "swc_refinement/final-world", "final-world")
    if voxel_dir is None or world_dir is None:
        logger.warning("Specimen-space inputs not found; skipping specimen resampling")
        return {"specimen_resampled": 0, "voxel_scale_um": None, "voxel_scale_source": None}

    # The transform stage carries acquisition.json forward; prefer it over inference.
    scale = scale_from_acquisition(DATA_DIR / "alignment" / "acquisition.json")
    source = "acquisition"
    if scale is None:
        scale = derive_scale(*matched_pair(voxel_dir, world_dir))
        source = "derived"
    staging = SCRATCH_DIR / "specimen_resampled_world"
    resample(world_dir, staging, spacing_um)

    inverse = tuple(1.0 / factor for factor in scale)
    written = 0
    for resampled in sorted(staging.glob("*.swc")):
        scale_swc(resampled, output_dir / resampled.name, inverse)
        written += 1
    shutil.rmtree(staging, ignore_errors=True)
    logger.info("Resampled %d specimen-space reconstruction(s) at %s um", written, spacing_um)
    return {
        "specimen_resampled": written,
        "voxel_scale_um": list(scale),
        "voxel_scale_source": source,
    }


def id_string(swc_path: Path) -> str:
    """Build the MouseLight ``idString`` for a reconstruction.

    Keeps the established ``<neuron>-<subject>`` form that downstream consumers expect,
    but parses the stem properly rather than taking the first two hyphen tokens.

    Parameters
    ----------
    swc_path : Path
        The reconstruction.

    Returns
    -------
    str
        The identifier, falling back to the bare stem if it does not parse.
    """
    try:
        parsed = parse_stem(swc_path.stem)
    except ReconstructionNameError:
        return swc_path.stem
    return f"{parsed.neuron_id}-{parsed.subject_id}"


def set_id_string(payload: object, value: str) -> int:
    """Rewrite every ``idString`` in a MouseLight payload.

    Parameters
    ----------
    payload : object
        The decoded JSON, at any depth.
    value : str
        The identifier to set.

    Returns
    -------
    int
        How many fields were rewritten.
    """
    count = 0
    if isinstance(payload, dict):
        for key, item in payload.items():
            if key == "idString":
                payload[key] = value
                count += 1
            else:
                count += set_id_string(item, value)
    elif isinstance(payload, list):
        for item in payload:
            count += set_id_string(item, value)
    return count


def annotate(swc_path: Path, destination: Path, mapper: CCFMorphologyMapper) -> None:
    """Annotate a CCF-space reconstruction and write it as MouseLight JSON.

    Parameters
    ----------
    swc_path : Path
        CCF-space reconstruction.
    destination : Path
        JSON file to write.
    mapper : CCFMorphologyMapper
        Maps coordinates to CCF structure ids.
    """
    morphology = read_swc(str(swc_path))
    mapper.annotate_morphology(morphology)
    destination.parent.mkdir(parents=True, exist_ok=True)
    MouseLightJsonWriter(morphology).write(str(destination))
    payload = json.loads(destination.read_text(encoding="utf-8"))
    value = id_string(swc_path)
    if set_id_string(payload, value) == 0 and isinstance(payload, dict):
        payload["idString"] = value
    destination.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def run() -> int:
    """Resample both coordinate spaces and annotate the CCF outputs.

    Returns
    -------
    int
        Process exit status. Non-zero when any reconstruction failed to annotate.
    """
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    args = parse_args()
    started = datetime.now(timezone.utc)

    aligned_dir = find_dir("alignment/aligned_swcs", "aligned_swcs")
    if aligned_dir is None:
        logger.error("No CCF-space reconstructions found under %s", DATA_DIR)
        return 1

    ccf_out = RESULTS_DIR / "final/ccf_space_reconstructions"
    resampled = SCRATCH_DIR / "aligned_resampled"
    resample(aligned_dir, resampled, args.spacing_um)

    mapper = CCFMorphologyMapper(
        reference_space_key="annotation/ccf_2017",
        resolution=int(args.ccf_resolution_um),
        cache_dir=str(SCRATCH_DIR),
    )
    failures: list[str] = []
    annotated = 0
    for swc_path in sorted(resampled.glob("*.swc")):
        destination = ccf_out / "swcs" / swc_path.name
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(swc_path, destination)
        try:
            annotate(destination, ccf_out / "jsons" / f"{swc_path.stem}.json", mapper)
            annotated += 1
        except Exception as error:  # noqa: BLE001 - one bad cell must not lose the run
            if args.fail_fast:
                raise
            logger.warning("Failed to annotate %s: %s", swc_path.name, error)
            failures.append(swc_path.stem)

    specimen = resample_specimen_space(
        args.spacing_um, RESULTS_DIR / "refinement/final-voxel-resampled"
    )

    write_stage_process(
        build_stage_process(
            STEP_NAME,
            resolve_code(
                "exaspim-swc-resample",
                url="https://github.com/peter-grotz/exaspim-swc-resample-capsule",
            ),
            start_time=started,
            output_path="final",
            parameters={
                "spacing_um": args.spacing_um,
                "ccf_resolution_um": args.ccf_resolution_um,
            },
            output_parameters={
                "ccf_swc_count": annotated + len(failures),
                "ccf_json_count": annotated,
                "failed": failures,
                **specimen,
            },
            experimenters=[e.strip() for e in args.experimenters.split(",") if e.strip()],
            notes=f"Resampled both coordinate spaces at {args.spacing_um} um.",
        ),
        RESULTS_DIR / "final",
    )

    logger.info("Annotated %d, failed %d", annotated, len(failures))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(run())
