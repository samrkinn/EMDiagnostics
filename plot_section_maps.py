import os
import sys
import numpy as np
import time
import renderapi
import json
import itertools as it
from argschema.fields import InputFile, InputDir, Str, Float, Int
from functools import partial
import pathos.multiprocessing as mp
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon
import matplotlib.cm as cm
from mpldatacursor import datacursor
import math

import mpld3
from mpld3 import plugins, utils


example = {
    "render":{
        "host":"http://em-131fs",
        "port":8080,
        "owner":"gayathri",
        "project":"MM2",
        "client_scripts":"/allen/programs/celltypes/workgroups/em-connectomics/gayathrim/nc-em2/Janelia_Pipeline/render_20170613/render-ws-java-client/src/main/scripts"
    },
    "input_stack":"mm2_acquire_8bit_reimage_Montage",
    "matchCollectionOwner":"gayathri_MM2",
    "matchCollection":"mm2_acquire_8bit_reimage_montage",
    "outdir":"/allen/programs/celltypes/workgroups/em-connectomics/gayathrim/nc-em2/Janelia_Pipeline/scratch/figures/montage",
    "zstart":1028,
    "zend":1028,
    "pool_size":20
}


def plot_section_maps(transformed_positions, tileIds, z, outdir, data=None):
    fig,ax = plt.subplots(1)

    #cmap = plt.get_cmap('RdYlBu')
    #nfloors = np.random.rand(len(transformed_positions))
    #colors = cmap(nfloors)
    print z
    cmap = plt.get_cmap('RdYlBu')
    #norm = mpl.colors.Normalize(vmin=min(data), vmax=max(data))
    norm = mpl.colors.Normalize(vmin=0, vmax=2)
    if data is None:
        colors = cmap(np.ones(len(transformed_positions)))
    else:
        colors = cmap(norm(data))

    patches = []

    for tp, ids in zip(transformed_positions, tileIds):
        polygon = Polygon(tp, closed=True, ec=(0,0,0,1), lw=2, label='${}$'.format(ids))
        patches.append(polygon)

    collection = PatchCollection(patches, alpha=0.5)
    cc = ax.add_collection(collection)
    collection.set_color(np.array(colors))
    ax.autoscale_view()
    plt.title(str(z))
    filename = os.path.join(outdir, '%s.png'%(str(z)))
    fig.savefig(filename, dpi=600)
    #plt.show()


def calculate_section_maps(render, input_stack, z):
    # read the tilespecs for the z
    tilespecs = render.run(renderapi.tilespec.get_tile_specs_from_z, input_stack, z)

    # for each tile specs compute the tile top left and bottom right corners
    transformed_positions = []
    for ts in tilespecs:
        x = 0
        y = 0
        #srcpts = [[x,y], [x+ts.width,y], [x+ts.width, y+ts.height], [x,y+ts.height]]
        srcpts = [[y,x], [y,x+ts.height], [y+ts.width, x+ts.height], [y+ts.width,x]]
        srcpts = np.array(srcpts)
        tforms = ts.tforms
        m = tforms.pop(0) # remove the lens correction transform

        dstpts = renderapi.transform.estimate_dstpts(tforms, src=srcpts)
        transformed_positions.append(dstpts)

    return transformed_positions

def get_tile_ids(tilespecs):
    tile_ids = [t.tileId for t in tilespecs]
    return tile_ids

def generate_section_maps(render, stack, outdir, z):
    transformed_positions = calculate_section_maps(render, stack, z)
    tilespecs = render.run(renderapi.tilespec.get_tile_specs_from_z, stack, z)
    tileIds = get_tile_ids(tilespecs)
    plot_section_maps(transformed_positions, tileIds, z, outdir)

if __name__ == "__main__":
    render = renderapi.Render(**example['render'])

    t1 = time.time()
    #plot_section_maps(transformed_positions, tileIds, z)

    mypartial = partial(generate_section_maps, render, example['input_stack'], example['outdir'])
    zvalues = range(example['zstart'], example['zend']+1)

    with mp.ProcessingPool(example['pool_size']) as pool:
        pool.map(mypartial, zvalues)
