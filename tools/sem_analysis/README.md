# SEM analysis Galaxy tool

This wrapper runs the locally available `sem-diameterj:1.0` Docker image. Copy
or link `sem_analysis.xml` into a Galaxy tool panel configuration and ensure the
image is available to every Galaxy job runner that can execute this tool.

Run the functional test with Planemo from this directory:

```bash
planemo test --docker sem_analysis.xml
```

The bundled test exercises the segmentation and Python QC path. The full mode
also launches DiameterJ 1.018 through the image's Xvfb-backed entrypoint.
