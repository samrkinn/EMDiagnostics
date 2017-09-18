
import renderapi
import os
import json
import numpy as np
from functools import partial

example = {
    "render":{
        "host": "http://em-131fs",
        "port": 8080,
        "owner": "gayathri",
        "project": "MM2",
        "client_scripts": "/allen/programs/celltypes/workgroups/em-connectomics/gayathrim/nc-em2/Janelia_Pipeline/render_20170613/render-ws-java-client/src/main/scripts"
    },
    "prestitched_stack": "mm2_acquire_8bit_reimage",
    "poststitched_stack": "mm2_acquire_8bit_reimage_Montage",
    "matchCollectionOwner": "gayathri_MM2",
    "matchCollection": "mm2_acquire_8bit_reimage_montage",
    "zstart": 1015,
    "zend": 1015,
    "pool_size": 20
}

#class DetectDisconnectedTilesParameters(RenderParameters):

def process_section(render, prestitched_stack, poststitched_stack, z):

    # get the tilespecs for both prestitched_stack and poststitched_stack
    pre_tilespecs = render.run(
                        renderapi.tilespec.get_tile_specs_from_z,
                        prestitched_stack,
                        z)
    post_tilespecs = render.run(
                        renderapi.tilespec.get_tile_specs_from_z,
                        poststitched_stack,
                        z)

    # insert all tiles from pre_tilespecs to Rtree
    [ridx.insert(i, (ts.minX, ts.minY, ts.maxX, ts.maxY)) for i, ts in enumerate(pre_tilespecs)]

    # pre tile_ids
    pre_tileIds = []
    [pre_tileIds.append(ts.tileId) for ts in pre_tilespecs]

    # post tile_ids
    post_tileIds = []
    [post_tileIds.append(ts.tileId) for ts in post_tileIds]

    missing_tileIds = list(set(pre_tileIds) - set(post_tileIds))
    print missing_tileIds
