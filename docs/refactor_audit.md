# Refactor audit

The original project archive contained exploratory notebooks, generated figures, compiled files, cache directories, checkpoint notebooks, external PDFs, and multiple exploratory model versions. This repository separates the reusable model implementation from legacy material.

Excluded from the public refactor by default:

- `.ipynb_checkpoints/`
- `__pycache__/`
- compiled files such as `.o`, `.pyc`, executable binaries
- external literature PDFs
- office lock files such as `.~lock.*`

Selected notebooks and original Python files are preserved under `legacy/selected_original_material/` for provenance.
