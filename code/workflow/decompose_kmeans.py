import numpy as np
import os
import sys
import pandas as pd
from datetime import datetime
from sklearn.cluster import KMeans
from pyminc.volumes.factory import *

# %% Script arguments
outdir = sys.argv[1]
decompdir = sys.argv[2]
matfile = sys.argv[3]
n_clusters = int(sys.argv[4])

# Testing arguments
# outdir = "data/outputs/isocortex_thalamus_left"
# decompdir = "data/outputs/isocortex_thalamus_left/decomposition_z/kmeans/10_clusters"
# matfile = "correlation_z.npy"
# n_clusters = 10

# %% Script specific calls

# KMeans: each elkan iteration in serial takes ~ 5 seconds
# 15 GB base + 15 GB per job
# Therefore running 5 jobs in parallel takes 90 GB, ~ 2 hours, for n_init=10, max_iter=500

model = KMeans(n_clusters=n_clusters, n_init=10, max_iter=500, verbose=1, n_jobs=5, algorithm="elkan", precompute_distances=False)
rscript = "code/workflow/decompose_kmeans.R"

# %% Create output directory
print("[{time}] Making output directory".format(time=datetime.now()))

os.makedirs(decompdir, exist_ok=True)


# %% Read data
print("[{time}] Reading data".format(time=datetime.now()))

mask_file_thalamus = os.path.join(outdir, "thalamus_mask_highres.mnc")
mask_file_cortex = os.path.join(outdir, "cortex_mask_lowres.mnc")
zmat_file = os.path.join(outdir, matfile)

mask_vol_thalamus = volumeFromFile(mask_file_thalamus, labels=True)
mask_vol_cortex = volumeFromFile(mask_file_cortex, labels=True)

p_thal = len(mask_vol_thalamus.data[mask_vol_thalamus.data > 0.5])
p_cortex = len(mask_vol_cortex.data[mask_vol_cortex.data > 0.5])

zmat = np.load(zmat_file)
zmat.shape

# %% Cluster/decompose
print("[{time}] Clustering/decomposing".format(time=datetime.now()))

model.fit(zmat)

print("[{time}] Done main computation".format(time=datetime.now()))

# %% Write out results
print("[{time}] Writing out npy arrays".format(time=datetime.now()))

np.save("{decompdir}/cluster_labels.npy".format(decompdir=decompdir), (model.labels_.astype('float64') + 1))
np.save("{decompdir}/cluster_centers.npy".format(decompdir=decompdir), model.cluster_centers_)

df = pd.DataFrame(data=model.get_params(), index=np.arange(0, 1))
df['inertia'] = model.inertia_
df.to_csv("{decompdir}/model_params.csv".format(decompdir=decompdir))

# %% Call Rscript to write out MNC files
print("[{time}] Writing out minc files via Rscript".format(time=datetime.now()))

cmd = "/usr/bin/time -v Rscript {rscript} {outdir} {decompdir} {n_clusters}".format(rscript=rscript, outdir=outdir, decompdir=decompdir, n_clusters=n_clusters)
os.system(cmd)

# %% Finally
print("[{time}] Done!".format(time=datetime.now()))
