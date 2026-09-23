"""Resample exaSPIM reconstructions and annotate the CCF-space ones.

Two resampling passes, both at 10 um:

* **CCF space** — the aligned reconstructions.
* **Specimen space** — resampled in *physical* coordinates and converted back to voxels.
  Specimen coordinates are voxels on an anisotropic grid, so resampling them directly
  would space nodes differently along each axis; see :mod:`scale`.

Code Ocean glue around ``neuron-tracing-utils``, whose resampling reads each SWC into an
SNT ``Tree`` and so needs a JVM and the Fiji jars. That stack keeps it out of the shared
library, which supplies the stage metadata record.
"""

import argparse
import logging
import os
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

from exaspim_swc_processing.resources import log_peak_memory
from exaspim_swc_processing.stage import (
    UPSTREAM_STAGES,
    build_stage_process,
    carry_forward,
    resolve_code,
    write_stage_process,
)
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
    parser.add_argument("--experimenters", default=os.environ.get("EXPERIMENTERS", ""))
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

    carried = carry_forward(DATA_DIR, RESULTS_DIR, UPSTREAM_STAGES)
    logger.info("Carried forward: %s", ", ".join(carried) or "nothing")

    aligned_dir = find_dir("alignment/aligned_swcs", "aligned_swcs")
    if aligned_dir is None:
        logger.error("No CCF-space reconstructions found under %s", DATA_DIR)
        return 1

    ccf_out = RESULTS_DIR / "final/ccf_space_reconstructions/swcs"
    resampled = SCRATCH_DIR / "aligned_resampled"
    resample(aligned_dir, resampled, args.spacing_um)

    ccf_out.mkdir(parents=True, exist_ok=True)
    ccf_count = 0
    for swc_path in sorted(resampled.glob("*.swc")):
        shutil.copy(swc_path, ccf_out / swc_path.name)
        ccf_count += 1

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
            parameters={"spacing_um": args.spacing_um},
            output_parameters={
                "ccf_swc_count": ccf_count,
                "stages_carried_forward": carried,
                **specimen,
            },
            experimenters=[e.strip() for e in args.experimenters.split(",") if e.strip()],
            notes=f"Resampled both coordinate spaces at {args.spacing_um} um.",
        ),
        RESULTS_DIR / "final",
    )

    logger.info("Resampled %d CCF-space reconstruction(s)", ccf_count)
    log_peak_memory()
    return 0


if __name__ == "__main__":
    sys.exit(run())
