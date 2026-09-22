# exaspim-swc-resample-capsule

Resamples exaSPIM reconstructions and annotates the CCF-space ones with CCF structure ids.

## Usage

Two resampling passes, both at 10 µm by default:

- **CCF space** — the aligned reconstructions, written to
  `/results/final/ccf_space_reconstructions/{swcs,jsons}` with MouseLight JSON alongside.
- **Specimen space** — written to `/results/refinement/final-voxel-resampled`.

Specimen coordinates are voxels on an anisotropic grid (`[0.748, 0.748, 1.0]` µm for
exaSPIM_794492), so resampling them directly would space nodes differently along each
axis. The pass therefore runs on the physical-coordinate reconstructions and converts the
result back to voxels. The scale is read from `acquisition.json`, which the transform stage needs for the
registration and carries forward. If it is absent, the same numbers are derived from a
matched voxel/world pair; that fallback is self-checking, since the same ratio must hold
on every axis for every node. Both paths were verified to agree on exaSPIM_794492.

A reconstruction that fails to annotate is logged and skipped; the failed stems are
recorded in the stage metadata and the capsule exits non-zero. Pass `--fail-fast` to abort
on the first failure.

## Level of Support

![support](https://img.shields.io/badge/support-supported-brightgreen)

## Installation

Stage metadata and stem parsing come from
[exaspim-swc-processing](https://github.com/peter-grotz/exaspim-swc-processing). The JVM
and pinned numeric stack stay here.
