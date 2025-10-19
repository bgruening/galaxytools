<h1>3Dtrees: Tile & Merge - Point cloud pre-processing for Model</h1>

**What it does**

This tool processes 3D point cloud data for tree segmentation by either:
- Tiling: Subsampling the input point cloud and creating tiles for processing
- Merging: Merging processed tiles back into the original point cloud resolution

**Tile Parameters**

- Tile Size: Size of tiles in meters (default: 50)
- Overlap: Overlap between tiles in meters (default: 20)
- Tiling Threshold: File size threshold in GB - files larger than this will be tiled (default: 3 GB)
- Points Threshold: Minimum number of points per tile - tiles with fewer points are deleted (default: 1000)
- Subsampling Resolution: Voxel size for subsampling in centimeters (default: 10 cm)

**Merge Parameters**

- Buffer Distance: Buffer distance for whole-tree assignment in meters (default: 0.2m)
- Minimum Cluster Size: Minimum points for a valid cluster (default: 300)
- Initial Search Radius: Starting radius for point reassignment in meters (default: 1.0m)
- Maximum Search Radius: Maximum radius for point reassignment in meters (default: 5.0m)
- Radius Step: Radius increment for expanding search in meters (default: 1.0m)
