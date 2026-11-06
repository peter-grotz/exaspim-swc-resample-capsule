# exaspim-swc-resample-capsule

Resamples exaSPIM reconstructions and annotates the CCF-space ones with CCF structure ids.

## Usage

Two resampling passes, both at 10 µm by default:

- **CCF space** — the aligned reconstructions, written to
  `/results/final/ccf_space_reconstructions/{swcs,jsons}` with MouseLight JSON alongside.
- **Specimen space** — the refined reconstructions in microns (`refinement/final-world`),
  written to `/results/refinement/final-world-resampled`.

A reconstruction that fails to annotate is logged and skipped; the failed stems are
recorded in the stage metadata and the capsule exits non-zero. Pass `--fail-fast` to abort
on the first failure.

## Level of Support

![support](https://img.shields.io/badge/support-supported-brightgreen)

## Installation

Stage metadata and stem parsing come from
[exaspim-swc-processing](https://github.com/AllenNeuralDynamics/exaspim-swc-processing). The JVM
and pinned numeric stack stay here.
