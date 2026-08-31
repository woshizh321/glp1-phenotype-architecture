# Reproducibility

1. Obtain official sources described in `DATA_SOURCES.md`.
2. Copy `config/paths.example.yaml` to `config/paths.yaml`.
3. Install packages in `environment/`.
4. Run `python run_pipeline.py --config config/paths.yaml --validate-inputs`.
5. Run `python run_pipeline.py --config config/paths.yaml`.

Generated row-level intermediates stay under the configured output root and must not be committed.
