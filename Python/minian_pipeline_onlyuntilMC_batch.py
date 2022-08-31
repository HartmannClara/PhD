# -*- coding: utf-8 -*-
"""
Created on Wed Aug 31 11:18:40 2022

@author: cha206

minian pipeline until saving of motion corrected video for further analysis in imagej
videocodec changed in visualization.py file to rawvideo so readable in imagej
"""
# %% loading modules
import itertools as itt
import os
import sys
import shutil

import holoviews as hv
import numpy as np
import xarray as xr
from dask.distributed import Client, LocalCluster
from holoviews.operation.datashader import datashade, regrid
from holoviews.util import Dynamic
#from IPython.core.display import display

#%% setting parameters

# Set up Initial Basic Parameters#
minian_path = "C:/Users/cha206/Anaconda3/envs/minian" 
#%% data loop


folderpath = "C:/Users/cha206/Data/DM/g2/2022_07_07/"

daylist = os.listdir(folderpath)
for x in daylist:
    dpath = folderpath + x +'/My_V4_Miniscope'
    minian_ds_path = os.path.join(dpath, "/minian")
    intpath = dpath + '/minian_intermediate'
        
    #dpath = "C:/Users/cha206/Data/DM/g2/2022_07_07/11_40_46/My_V4_Miniscope"
    #minian_ds_path = os.path.join(dpath, "minian")
    #intpath = "C:/Users/cha206/Data/DM/g2/2022_07_07/11_40_46/My_V4_Miniscope/minian_intermediate"
    subset = dict(frame=slice(0, None))
    subset_mc = None
    interactive = True
    output_size = 100
    n_workers = int(os.getenv("MINIAN_NWORKERS", 4))
    
    param_save_minian = {
        "dpath": minian_ds_path,
        "meta_dict": dict(session=-2, animal=-5),#change so metadata is saved correctly!!!!
        "overwrite": True,
    }
    
    #Pre-processing Parameters
    param_load_videos = {
        "pattern": '0_corr_crop+\.avi$'        , #"[0-0]+\.avi$",
        "dtype": np.uint8,
        "downsample": dict(frame=1, height=1, width=1),
        "downsample_strategy": "subset",
    }
    param_denoise = {"method": "median", "ksize": 7}
    param_background_removal = {"method": "tophat", "wnd": 15}
    
    # Motion Correction Parameters#
    subset_mc = None
    param_estimate_motion = {"dim": "frame"}
    
    # Initialization Parameters#
    # =============================================================================
    # param_seeds_init = {
    #     "wnd_size": 1000,
    #     "method": "rolling",
    #     "stp_size": 500,
    #     "max_wnd": 15,
    #     "diff_thres": 3,
    # }
    # param_pnr_refine = {"noise_freq": 0.06, "thres": 1}
    # param_ks_refine = {"sig": 0.05}
    # param_seeds_merge = {"thres_dist": 10, "thres_corr": 0.8, "noise_freq": 0.06}
    # param_initialize = {"thres_corr": 0.8, "wnd": 10, "noise_freq": 0.06}
    # param_init_merge = {"thres_corr": 0.8}
    # 
    # # CNMF Parameters#
    # param_get_noise = {"noise_range": (0.06, 0.5)}
    # param_first_spatial = {
    #     "dl_wnd": 5,
    #     "sparse_penal": 0.01,
    #     "update_background": True,
    #     "size_thres": (25, None),
    # }
    # param_first_temporal = {
    #     "noise_freq": 0.06,
    #     "sparse_penal": 1,
    #     "p": 1,
    #     "add_lag": 20,
    #     "jac_thres": 0.2,
    # }
    # param_first_merge = {"thres_corr": 0.8}
    # param_second_spatial = {
    #     "dl_wnd": 5,
    #     "sparse_penal": 0.01,
    #     "update_background": True,
    #     "size_thres": (25, None),
    # }
    # param_second_temporal = {
    #     "noise_freq": 0.06,
    #     "sparse_penal": 1,
    #     "p": 1,
    #     "add_lag": 20,
    #     "jac_thres": 0.4,
    # }
    # 
    # =============================================================================
    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["MINIAN_INTERMEDIATE"] = intpath
    
    
    # %% import minian
    
    sys.path.append(minian_path)
    from minian.cnmf import (
        compute_AtC,
        compute_trace,
        get_noise_fft,
        smooth_sig,
        unit_merge,
        update_spatial,
        update_temporal,
    )
    from minian.initialization import (
        gmm_refine,
        initA,
        initbf,
        initC,
        intensity_refine,
        ks_refine,
        pnr_refine,
        seeds_init,
        seeds_merge,
    )
    from minian.motion_correction import apply_transform, estimate_motion
    from minian.preprocessing import denoise, remove_background
    from minian.utilities import (
        TaskAnnotation,
        get_optimal_chk,
        load_videos,
        open_minian,
        save_minian,
    )
    from minian.visualization import (
        CNMFViewer,
        VArrayViewer,
        generate_videos,
        visualize_gmm_fit,
        visualize_motion,
        visualize_preprocess,
        visualize_seeds,
        visualize_spatial_update,
        visualize_temporal_update,
        write_video,
    )
    #%% module initialization
    dpath = os.path.abspath(dpath)
    #hv.notebook_extension("bokeh", width=100) #oly for jup notebook?
    
    #%% start cluster, not needed
    
    # =============================================================================
    # cluster = LocalCluster(
    #     n_workers=n_workers,
    #     memory_limit="3GB",#2 is minimum
    #     resources={"MEM": 1},
    #     threads_per_worker=2,
    #     dashboard_address=":8787",
    # )
    # annt_plugin = TaskAnnotation()
    # cluster.scheduler.add_plugin(annt_plugin)
    # client = Client(cluster)
    # =============================================================================
    
    
    #%%loading and chunking
    
    varr = load_videos(dpath, **param_load_videos)
    chk, _ = get_optimal_chk(varr, dtype=float)
    
    #save as array
    varr = save_minian(
        varr.chunk({"frame": chk["frame"], "height": -1, "width": -1}).rename("varr"),
        intpath,
        overwrite=True,
    )
    
    
    varr_ref = varr.sel(subset)#nothing for analysis subsetting
    
    #%% glow removal, denoising and background removal
    
    varr_min = varr_ref.min("frame").compute()
    varr_ref = varr_ref - varr_min
    
    varr_ref = denoise(varr_ref, **param_denoise)
    
    varr_ref = remove_background(varr_ref, **param_background_removal)
    
    varr_ref = save_minian(varr_ref.rename("varr_ref"), dpath=intpath, overwrite=True)
    #%% motion correction
    motion = estimate_motion(varr_ref.sel(subset_mc), **param_estimate_motion)
    #motion = save_minian(
    #    motion.rename("motion").chunk({"frame": chk["frame"]}), **param_save_minian)
    
    Y = apply_transform(varr_ref, motion, fill=0)
    
    Y_fm_chk = save_minian(Y.astype(float).rename("Y_fm_chk"), intpath, overwrite=True)
    Y_hw_chk = save_minian(
        Y_fm_chk.rename("Y_hw_chk"),
        intpath,
        overwrite=True,
        chunks={"frame": -1, "height": chk["height"], "width": chk["width"]},
    )
    
    
    #%% save mc video
    write_video(Y_fm_chk, "minian_mc.avi", dpath)
    



