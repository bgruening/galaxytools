The fixtures contain 500 points and SAT/ForestMamba prediction columns for item
11353, generation `0b0a6f1e-6a4e-4d8e-9f3c-1d2e3f4a5b6c`.

- `11353` and `base.laz` are identical standardized inputs. The extensionless
  name supplies the item identity to Galaxy's dataset wrapper.
- `sat/11353` and `sat.laz` contain the SAT predictions.
- `legacy/11353` contains SAT and ForestMamba predictions together.
- `sat.parquet` and `forestmamba.parquet` are the corresponding indexed
  sidecars, including item and generation metadata.
- `member_both.laz` is the assembled base with both methods. The golden files
  match the published 0.1.0 container built from upstream commit
  `28eb4811ef34901619a9a26c406a074cdb5fa7ea`.
- `11354` duplicates `base.laz` under a different item identity, to test
  rejection of mismatched sidecars.
- `sat_short.parquet` has the first 499 rows of `sat.parquet` and a manifest
  point count of 499, to test point-count validation.
- `member_changed.laz` duplicates `member_both.laz` with
  `PredInstance_SAT[123]` incremented by one, to test value-level checking.

The tests compare binary products exactly and inspect parsed JSON reports.
Failure tests require exit code 2 and the expected rejection message, so an
unrelated crash cannot pass as a successful refusal.
