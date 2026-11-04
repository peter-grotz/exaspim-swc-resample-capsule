"""Resample exaSPIM reconstructions and annotate the CCF-space ones.

Two resampling passes, both at 10 um:

* **CCF space** — the aligned reconstructions.
* **Specimen space** — the refined reconstructions in physical coordinates.

"""

import argparse
import logging
import os
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

from exaspim_swc_processing.stage import (
    UPSTREAM_STAGES,
    build_stage_process,
    carry_forward,
    resolve_code,
    write_stage_process,
)

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


def resample_specimen_space(spacing_um: float, output_dir: Path) -> int:
    """Resample the physical-coordinate specimen-space reconstructions.

    Parameters
    ----------
    spacing_um : float
        Target spacing between nodes, in microns.
    output_dir : Path
        Directory to write the resampled reconstructions to.

    Returns
    -------
    int
        How many were written.
    """
    world_dir = find_dir("refinement/final-world", "swc_refinement/final-world", "final-world")
    if world_dir is None:
        logger.warning("Specimen-space inputs not found; skipping specimen resampling")
        return 0
    resample(world_dir, output_dir, spacing_um)
    written = len(list(output_dir.glob("*.swc")))
    logger.info("Resampled %d specimen-space reconstruction(s) at %s um", written, spacing_um)
    return written


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

    resample_specimen_space(args.spacing_um, RESULTS_DIR / "refinement/final-world-resampled")
    write_stage_process(
        build_stage_process(
            STEP_NAME,
            resolve_code(url="https://github.com/peter-grotz/exaspim-swc-resample-capsule"),
            started,
            experimenters=[e.strip() for e in args.experimenters.split(",") if e.strip()],
        ),
        RESULTS_DIR / "final",
    )

    logger.info("Resampled %d CCF-space reconstruction(s)", ccf_count)
    return 0


if __name__ == "__main__":
    sys.exit(run())
