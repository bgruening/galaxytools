# RayCloudTools comparison fixtures

`mikro` is the no-GPS LAZ input.
The file is deliberately extensionless, matching the original wrapper-test
convention: the unchanged command appends `.laz` to the Galaxy dataset name.

Expected outputs were generated from the 12,917-point fixture with the released
`ghcr.io/3dtrees-earth/3dtrees_rct:1.2.2` image, digest
`sha256:ebc6e35dcd62af4a43272f07ed4c82a7281a46bc33705ed23c1834fa3f858268`.
Tests use the wrapper defaults, then Treeinfo, then Treeinfo plus branch
segmentation and branch data. The branch case has separate LAZ and TXT references.
Both PLY references are shared across all three cases.

TXT outputs use content comparison. LAZ/PLY outputs use size comparison, as in
the previous tests, alongside explicit label, ownership and mesh-header checks;
size comparison is not a byte-for-byte geometry or point-value comparison.

Known producer limitation: dotted dataset identifiers can produce an incorrect
Treeinfo path in v1.2.2. These fixtures preserve the legacy supported naming
convention; this test-only change does not repair that runtime limitation.
