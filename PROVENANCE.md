# Provenance and publication scope

This repository is the public companion to Youngchae Lee's butane computational chemistry coursework extension. Publication was requested by the owner on 2026-10-15. Code, input preparation, log auditing and documentation were developed with AI assistance. Numerical claims here come from actually executed calculations and the supplied analysis, with approximations stated explicitly.

## Data lineage

- `data/original_coursework.csv`: four original user coursework logs dated 2025-09-28. The immutable raw logs are retained in the private project archive. No values were recovered by fitting the printed report.
- `data/minima_audit.csv`: three newly executed Gaussian 16 C.01 Opt/Freq jobs on 2026-10-15. Final electronic energy and geometry belong to the completed optimization; printed thermal/Gibbs sums belong to the subsequent frequency calculation. Gibbs sums are rounded to six decimal places in Gaussian output, limiting energy differences to approximately 0.0026 kJ/mol precision from rounding alone.
- `data/torsion_audit.csv`: seven new Gaussian 16 C.01 constrained optimizations on 2026-10-15. All completed with exit code zero. Source log hashes were compared between compute host and local retrieval. Energy extraction uses the last SCF preceding the completed optimization and its corresponding orientation table.
- `data/input_provenance.csv`: input names, original geometry source hashes, rigid transformations, target angles and roles fixed before calculation. Rotations preserve each fragment's geometry about the central C5–C8 axis.
- `data/final_geometries.json`: optimized coordinates in angstrom from the seven completed torsion jobs. Atom ordering matches `inputs/constrained/`.
- `data/input_hashes.json`: hashes of the exact published new inputs, identical to those executed.
- `results/`: derived outputs from `analyze.py`. Curve values between calculations are predictions, not independent reference energies.

The fit/held-out angle assignment and model orders were recorded before running the seven jobs. Mirroring a geometry or duplicating a log is not an independent validation sample. The minimum calculations had already been inspected before this split was designed and are labeled diagnostics.

## Reproducibility boundary

The public package reproduces numerical analysis from audited extracts, and includes inputs for rerunning the new calculations. Raw source logs, historical input files and checkpoints are not distributed in this release. Consequently a reader cannot independently authenticate the entire historical extraction from this repository alone; hashes identify the retained evidence but do not replace it. Public final geometries and condition/convergence columns allow inspection of the extracted records.

Personal paths, hostnames, server accounts, credentials, license text, report/presentation PDFs, private search inventories and other students' reference material are excluded. The package contains no Gaussian executable or source code. No claim of independent human review of every generated line is made. No automated software test suite or hosted CI result is claimed in this release; scientific log/convergence/geometry audits and the analysis were actually executed.

## Continuing publication

Future updates should publish reviewed inputs, curated extracts, reproducible analysis and explicit limitations at this repository. Preserve earlier evaluated datasets and keep subsequent fit/validation splits identifiable. Do not upload the private working directory wholesale. Public Git history is the record of published releases, while complete raw backups remain in the private archive.
