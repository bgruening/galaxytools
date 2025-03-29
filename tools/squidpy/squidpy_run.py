import scanpy as sc
import anndata as ad
import squidpy as sq
import pandas as pd
import matplotlib.pyplot as plt
import argparse
from PIL import Image

parser = argparse.ArgumentParser(description="Command parameters for squidpy", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-a", "--adata", help="Anndata file")
parser.add_argument("-l", "--library_id", help="Visium library ID. Defaults to None")
args = parser.parse_args()

def readin(anndata, library_id):
    adata=ad.read_h5ad(anndata)
    scale = adata.uns['spatial'][library_id]['scalefactors']['tissue_hires_scalef']
    img = sq.im.ImageContainer(adata.uns['spatial'][library_id]['images']['hires'],
                            scale=scale, library_id=library_id)
    return adata, img

def cluster_features(features: pd.DataFrame, like=None):
    """
    Calculate leiden clustering of features.
    Specify filter of features using `like`.
    """
    # filter features
    if like is not None:
        features = features.filter(like=like)
    # create temporary adata to calculate the clustering
    adata = ad.AnnData(features)
    # important - feature values are not scaled, so need to scale them before PCA
    sc.pp.scale(adata)
    # calculate leiden clustering
    sc.pp.pca(adata, n_comps=min(10, features.shape[1] - 1))
    sc.pp.neighbors(adata)
    sc.tl.leiden(adata) 

    return adata.obs["leiden"]


img_count = 1
adata_file = args.adata
library_id = args.library_id

#import data
adata, img = readin(adata_file, library_id)

# # Visualize cluster annotation in spatial context
# cluster_image = sq.pl.spatial_scatter(adata, color="cluster")

# # Show all image channels 
# channel_image = img.show(channelwise=True)

# Image process by smoothing and segmentation
sq.im.process(img=img, layer="image", method="smooth")
sq.im.segment(img=img, layer="image_smooth", method="watershed", channel=0, ax=None, chunks=1000)

# plot the resulting segmentation
fig, ax = plt.subplots(1, 2)

#HAND TO USER FOR CROP

img.show(
    layer="image", 
    channel=0, 
    ax=None,
    save="img_" + str(img_count) + ".png",)
img_count+=1
img.show(
    layer="segmented_watershed",
    channel=0,
    ax=None,
   save="img_" + str(img_count) + ".png",)
img_count+=1


# define image layer to use for segmentation
features_kwargs = {"segmentation": {"label_layer": "segmented_watershed"}}
# calculate segmentation features
sq.im.calculate_image_features(
    adata,
    img,
    features="segmentation",
    layer="image",
    key_added="features_segmentation",
    n_jobs=None,
    features_kwargs=features_kwargs,
)
# plot results and compare with gene-space clustering
features = adata.obsm["features_segmentation"]
features = [i for i in features if "intensity_mean" in i]
features.append("segmentation_label")
features.append("clusters")


feature_subset = sq.pl.extract(adata, "features_segmentation")
for i in features:
    sq.pl.spatial_scatter(
        feature_subset,
        color=i,
        frameon=True,
        ncols=1,
        ax=None,
        save="img_" + str(img_count) + ".png",
    )
    img_count+=1



# # define different feature calculation combinations
params = {
    # all features, corresponding only to tissue underneath spot
    # "features_orig": {
    #     "features": ["summary", "texture", "histogram"],
    #     "scale": 1.0,
    #     "mask_circle": True,
    # },
    # summary and histogram features with a bit more context, original resolution
    "features_context": {"features": ["summary", "histogram"], "scale": 1.0},
    # summary and histogram features with more context and at lower resolution
    
    # "features_lowres": {"features": ["summary", "histogram"], "scale": 0.25},
}

for feature_name, cur_params in params.items():
    # features will be saved in `adata.obsm[feature_name]`
    sq.im.calculate_image_features(adata, img, layer="image", key_added=feature_name, n_jobs=1, **cur_params)

# combine features in one dataframe
adata.obsm["features"] = pd.concat([adata.obsm[f] for f in params.keys()], axis="columns")

# make sure that we have no duplicated feature names in the combined table
adata.obsm["features"].columns = ad.utils.make_index_unique(adata.obsm["features"].columns)


adata.obs["features_summary_cluster"] = cluster_features(adata.obsm["features"], like="summary")
adata.obs["features_histogram_cluster"] = cluster_features(adata.obsm["features"], like="histogram")
# adata.obs["features_texture_cluster"] = cluster_features(adata.obsm["features"], like="texture")

sc.set_figure_params(facecolor="white", figsize=(8, 8))
sq.pl.spatial_scatter(
    adata,
    color=[
        "features_summary_cluster",
        "features_histogram_cluster",
        # "features_texture_cluster",
        "clusters",
    ],
    ncols=3,
    save="img_" + str(img_count) + ".png",
)
img_count+=1

# # Neightborhood enrichment analsyis
sq.gr.spatial_neighbors(adata)
sq.gr.nhood_enrichment(adata, cluster_key="clusters")
# # plot neighborhood enrichment
sq.pl.nhood_enrichment(adata, cluster_key="clusters", save="img_" + str(img_count) + ".png")
img_count+=1

# # Co-occurrence across spaital dimensions
sq.gr.co_occurrence(adata, cluster_key="clusters")
sq.pl.co_occurrence(
    adata,
    cluster_key="clusters",
    clusters=['0'],
    figsize=(8, 4),
)


#plot for every single cluster
clusters_list = adata.obs['clusters'].cat.categories
for i in clusters_list:
    img_file = "img_" + str(img_count) + ".png"
    sq.gr.co_occurrence(adata, cluster_key="clusters")
    sq.pl.co_occurrence(
        adata,
        cluster_key="clusters",
        clusters=i,
        figsize=(8, 4),
        save=img_file,
    )
    img_count += 1

img_files = []
for i in range(1,img_count):
    current = "img_" + str(i) + ".png"
    img_files.append(current)

images = [
    Image.open("./figures/" + f)
    for f in img_files
    ]

path = "./combined.pdf"

new_images = []

for file in images:
    if file.mode == 'RGBA':
        file = file.convert('RGB')
    new_images.append(file)

new_images[0].save(path, "PDF", resolution=100.0, save_all=True, append_images=new_images[1:])


# # Ligand Receptor analysis 
# sq.gr.ligrec(
#     adata,
#     n_perms=100,
#     cluster_key="cluster",
# )
# sq.pl.ligrec(
#     adata,
#     cluster_key="cluster",
#     source_groups="Hippocampus",
#     target_groups=["Pyramidal_layer", "Pyramidal_layer_dentate_gyrus"],
#     means_range=(3, np.inf),
#     alpha=1e-4,
#     swap_axes=True,
# )

# # spatially variable genes 
# genes = adata[:, adata.var.highly_variable].var_names.values[:1000]
# sq.gr.spatial_autocorr(
#     adata,
#     mode="moran",
#     genes=genes,
#     n_perms=100,
#     n_jobs=1,
# )

# sq.pl.spatial_scatter(adata, color=["Olfm1", "Plp1", "Itpka", "cluster"])