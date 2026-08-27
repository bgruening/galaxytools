# SEM analysis Galaxy tool

This wrapper runs
`quay.io/galaxy/sem-analysis-fiji-diameterj:0.2`. Copy or link
`sem_analysis.xml` into a Galaxy tool panel configuration and ensure the image
is available to every Galaxy job runner that can execute this tool.

Run the functional test with Planemo from this directory:

```bash
planemo test --docker sem_analysis.xml
```

The bundled tests exercise the Fiji segmentation and Python-QC paths. Full mode
launches DiameterJ 1.018 through the image's Xvfb-backed entrypoint.
