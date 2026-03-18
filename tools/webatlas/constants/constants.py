#!/usr/bin/env python3
from vitessce import (
    DataType as dt,
    FileType as ft,
    Component as cm,
    CoordinationType as ct,
)

SINGLE_ZARR = "anndata.zarr"
OBS = "obs"

# @TODO: remove deprecated and unused data types
# Data types with ordered file types
DATA_TYPES = {
    OBS: [
        (SINGLE_ZARR, ft.ANNDATA_ZARR),
    ],
    # dt.OBS_EMBEDDING: [
    #     ("obsEmbedding.anndata.zarr", ft.OBS_EMBEDDING_ANNDATA_ZARR),
    #     ("obsEmbedding.csv", ft.OBS_EMBEDDING_CSV),
    # ],
    # dt.OBS_LABELS: [
    #     ("obsLabels.anndata.zarr", ft.OBS_LABELS_ANNDATA_ZARR),
    #     ("obsLabels.csv", ft.OBS_LABELS_CSV),
    # ],
    # dt.OBS_LOCATIONS: [
    #     ("obsLocations.anndata.zarr", ft.OBS_LOCATIONS_ANNDATA_ZARR),
    #     ("obsLocations.csv", ft.OBS_LOCATIONS_CSV),
    # ],
    dt.MOLECULES: [
        ("molecules.json", ft.MOLECULES_JSON),
    ],
    # dt.OBS_SETS: [
    #     ("cell-sets.json", ft.CELL_SETS_JSON),
    #     ("obsSets.anndata.zarr", ft.OBS_SETS_ANNDATA_ZARR),
    #     ("anndata-cell-sets.zarr", ft.ANNDATA_CELL_SETS_ZARR),
    #     (SINGLE_ZARR, ft.OBS_SETS_ANNDATA_ZARR),
    # ],
    dt.RASTER: [
        ("raster.ome-zarr", "raster.ome-zarr"),
        ("raster.json", ft.RASTER_JSON),
    ],
    # dt.OBS_FEATURE_MATRIX: [
    #     ("obsFeatureMatrix.anndata.zarr", ft.OBS_FEATURE_MATRIX_ANNDATA_ZARR),
    #     ("expression-matrix.zarr", ft.EXPRESSION_MATRIX_ZARR),
    #     ("anndata-expression-matrix.zarr", ft.ANNDATA_EXPRESSION_MATRIX_ZARR),
    #     ("clusters.json", ft.CLUSTERS_JSON),
    #     ("genes.json", ft.GENES_JSON),
    #     (SINGLE_ZARR, ft.OBS_FEATURE_MATRIX_ANNDATA_ZARR),
    # ],
    dt.NEIGHBORHOODS: [
        ("neighborhoods.json", ft.NEIGHBORHOODS_JSON),
    ],
    dt.GENOMIC_PROFILES: [
        ("genomic-profiles.zarr", ft.GENOMIC_PROFILES_ZARR),
    ],
}

DEFAULT_OPTIONS = {
    ft.ANNDATA_ZARR: {
        "mappings": {
            "obsm/X_umap": [0, 1],
        },
        "factors": [
            "obs/sample",
        ],
        "spatial": {
            "xy": "obsm/spatial",
        },
        "sets": ["obs/sample"],
        "matrix": "X",
    }
}


def hconcat(*cms):
    return "({})".format(("|").join(cms))


def vconcat(*cms):
    return "({})".format(("/").join(cms))


DEFAULT_LAYOUTS = {
    "minimal": hconcat(cm.SPATIAL.value, cm.LAYER_CONTROLLER.value),
    "simple": hconcat(
        cm.SPATIAL.value,
        hconcat(
            cm.LAYER_CONTROLLER.value, vconcat(cm.FEATURE_LIST.value, cm.OBS_SETS.value)
        ),
    ),
    "advanced": hconcat(
        cm.LAYER_CONTROLLER.value,
        cm.SPATIAL.value,
        hconcat(
            vconcat(cm.SCATTERPLOT.value, cm.OBS_SETS.value),
            cm.FEATURE_LIST.value,
        ),
        cm.GENOMIC_PROFILES.value,
    ),
}

# Coordination Types required by Components/Views
COMPONENTS_COORDINATION_TYPES = {cm.SCATTERPLOT: [ct.EMBEDDING_TYPE]}

# Data Types required by Components/Views
COMPONENTS_DATA_TYPES = {
    cm.SCATTERPLOT: set([dt.OBS_EMBEDDING]),
    cm.HEATMAP: set([dt.OBS_FEATURE_MATRIX]),
    cm.SPATIAL: set([dt.RASTER, dt.OBS_LOCATIONS, dt.MOLECULES]),
    cm.LAYER_CONTROLLER: set([dt.RASTER, dt.OBS_LOCATIONS, dt.MOLECULES]),
    cm.GENOMIC_PROFILES: set([dt.GENOMIC_PROFILES]),
    cm.FEATURE_LIST: set([dt.OBS_FEATURE_MATRIX]),
    cm.OBS_SETS: set([dt.OBS_SETS]),
    cm.OBS_SET_SIZES: set([dt.OBS_SETS]),
    cm.OBS_SET_FEATURE_VALUE_DISTRIBUTION: set([dt.OBS_FEATURE_MATRIX]),
}
