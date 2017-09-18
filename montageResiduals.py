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
    "zstart":1015,
    "zend":1015,
    "pool_size":20
}

def plot_section_maps(transformed_positions, data):
    fig,ax = plt.subplots(1)

    #cmap = plt.get_cmap('RdYlBu')
    #nfloors = np.random.rand(len(transformed_positions))
    #colors = cmap(nfloors)

    cmap = cm.cool
    #norm = mpl.colors.Normalize(vmin=min(data), vmax=max(data))
    norm = mpl.colors.Normalize(vmin=0, vmax=2)
    colors = cmap(norm(data))

    patches = []

    for tp in transformed_positions:
        polygon = Polygon(tp, closed=True, ec=(0,0,0,1), lw=2)
        patches.append(polygon)

    collection = PatchCollection(patches)
    ax.add_collection(collection)
    collection.set_color(np.array(colors))
    ax.autoscale_view()
    plt.show()



def plot_stack_statistics(data, zvalues):
    # scatter plot for plotting group_mean or group_variance or group_std_deviation
    # input - a list of numbers representing the data

    # Using a closure to access data. Ideally you'd use a "functor"-style class.
    def formatter(**kwargs):
        dist = abs(np.array(x) - kwargs['x'])
        i = dist.argmin()
        return '\n'.join(str(x[i]))

    if len(data) != len(zvalues):
        print "Missing data or zvalues for plotting"
    else:
        fig, ax = plt.subplots()
        x = np.array(zvalues)
        y = np.array(data)
        colors = np.random.rand(len(data))
        ax.scatter(x,y, c=colors, alpha=0.5)
        datacursor(hover=True, formatter=formatter)
        plt.show()



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

def get_section_ids(render, input_stack):
    data = render.run(renderapi.stack.get_stack_sectionData, input_stack)
    zvalues = []
    sectionIds = []
    for d in data:
        zvalues.append(d['z'])
        sectionIds.append(d['sectionId'])
    return zvalues, sectionIds

def get_groups(render, matchCollectionOwner, matchCollection):
    groups = render.run(renderapi.pointmatch.get_match_groupIds, matchCollection, owner=matchCollectionOwner)
    groups = np.array(map(int, groups))
    groups.sort()
    return groups



def get_tile_ids(tilespecs):
    tile_ids = [t.tileId for t in tilespecs]
    return tile_ids


def stats_and_outliers(data, ids, cutoff = 1, method = 'std', greater_than = True, normalizer = 1, **kwargs):
    '''Returns dict of various statistics and information about outliers given
    an input dataset. The input data in this case is usually a list of tile
    ratios (area and perimeter, for example) or list of tile residuals (from
    point matches)

    cutoff -- the value to compare data to determine whether or not it is
    an outlier. Meaning of this number will depend on 'method'. If using 'fixed'
    cutoff, this will assume that the input is a ratio. May need to subtract 1
    from cutoff to use it properly in that case. Can specify a positive and neg
    cutoff, if so cutoff will be a 2-element list with the first element as
    the positive cutoff
    method -- 'std' or 'fixed'. 'std' judges outliers based on how many standard
    deviations data is from the mean, where cutoff is the number of stdevs.
    'fixed' uses an absolute cutoff, specified by cutoff
    greater_than -- bool, if True, will only tag a tile as an outlier if it is
    greater than the cutoff (pos. stdev away from the mean, for example) but
    not lower than it.
    normalizer -- keep this at 1 most of the time, specifies whether variance
    and stdev are calculated using N-(normalizer) normalization
    '''

    output = {}

    mean = np.nanmean(data)
    var = np.nanvar(data, ddof = normalizer)
    med = np.nanmedian(data)
    std = var ** (0.5)

    if method == 'fixed':
        if greater_than:
            outlier_idx =  [i for i in range(len(data)) if data[i]-1 > cutoff]
        else:
            if type(cutoff) != list:
                outlier_idx = [i for i in range(len(data)) if abs(data[i] - 1) >= cutoff]
            else:
                outlier_idx = [i for i in range(len(data)) if data[i] - 1 >= cutoff[0] or data[i] - 1 <= cutoff[2]]
    elif method == 'std':
        if greater_than:
            outlier_idx = [i for i in range(len(data)) if data[i] - mean >= cutoff*std]
        else:
            outlier_idx = [i for i in range(len(data)) if abs(data[i] - mean) >= cutoff*std]

    # find ids of tiles corresponding to those indecies
    outlier_ids = ids[outlier_idx]
    outlier_count = len(outlier_idx)
    outlier_percent = float(outlier_count)/len(data) * 100.0

    output['mean'] = mean
    output['variance'] = var
    output['median'] = med
    output['stdev'] = std
    output['outlier ids'] = outlier_ids
    output['outlier idx'] = outlier_idx
    output['outlier count'] = outlier_count
    output['outlier percent'] = outlier_percent

    return output


def calculate_residual_for_match(match, p_tiles_id, p_tiles, q_tiles_id = None, q_tiles = None):
    '''Returns the residual for a specified point match.

    match -- an dict from a list of point matches specifying the corresponding
    points and the tiles of interest
    p_tiles_id -- list or array of tile ids contained in p_tiles
    p_tiles -- list of tilespec objects
    q_tiles_id/q_tiles -- same function as the above two, may be the same.
    Different if using cross-section point matches.
    '''

    if q_tiles_id == None:
        q_tiles_id = p_tiles_id
    if q_tiles == None:
        q_tiles = p_tiles

    p_points = np.column_stack((np.array(match['matches']['p'][0]), np.array(match['matches']['p'][1])))
    q_points = np.column_stack((np.array(match['matches']['q'][0]), np.array(match['matches']['q'][1])))

    p_idx = np.where(p_tiles_id == match['pId'])
    q_idx = np.where(q_tiles_id == match['qId'])

    # skip when match is out of stack? bug?
    # the render pm collections had some matches to tiles not found in the
    # stack tilespecs, but ignoring them made results match perfectly with
    # EM_aligner matlab code
    if p_idx[0].size == 0 and q_idx[0].size == 0:
        return 'flag', match['pId'], match['qId']
    elif q_idx[0].size == 0:
        return 'flag', [], match['qId']
    elif p_idx[0].size ==0:
        return 'flag', match['pId'], []

    p_idx = p_idx[0][0]
    q_idx = q_idx[0][0]
    # only works for PolynomialTransform2D or Affine classes, has to be second entry in tforms
    if p_tiles[p_idx].tforms[1].className == 'mpicbg.trakem2.transform.PolynomialTransform2D' and q_tiles[q_idx].tforms[1].className == 'mpicbg.trakem2.transform.PolynomialTransform2D':
        p_points_tformed = tform_poly2d_asAffine(p_points, p_tiles[p_idx].tforms[1].params)
        q_points_tformed = tform_poly2d_asAffine(q_points, q_tiles[q_idx].tforms[1].params)
    elif p_tiles[p_idx].tforms[1].className == 'mpicbg.trakem2.transform.AffineModel2D' and q_tiles[q_idx].tforms[1].className == 'mpicbg.trakem2.transform.AffineModel2D':
        p_points_tformed = p_tiles[p_idx].tforms[1].tform(p_points)
        q_points_tformed = q_tiles[q_idx].tforms[1].tform(q_points)
    else:
        print 'Could not compatible transform 2D in second transform slot, no transformation taken'

    sub = p_points_tformed - q_points_tformed

    # equivalent to mean of all sqrt((px_n - qy_n)^2 + (py_n - qy_n)^2) for all n
    residual = np.mean(np.sqrt(np.add(np.multiply(sub[:,0],sub[:,0]),np.multiply(sub[:,1],sub[:,1]))))
    return residual, p_idx, q_idx

def get_matches(render, stack, matchCollectionOwner, matchCollection, pgroupId, qgroupId):
    allmatches = render.run(
                    renderapi.pointmatch.get_macthes_from_group_to_group,
                    matchCollection,
                    pgroupId,
                    qgroupId,
                    owner=matchCollectionOwner)

    return allmatches

def area_of_polygon(x,y):
    # solution is based on shoelace formula
    return 0.5 * np.abs(np.dot(x, np.roll(y,1)) - np.dot(y, np.roll(x,1)))


def compute_mean_tile_residuals(tile_residuals):
    tile_mean = {}
    unconnected_tile_ids = []

    #print tile_residuals.keys()

    # loop over each tile and compute the mean residual for each tile
    # iteritems is specific to py2.7
    for key, residuals in tile_residuals.iteritems():
        if tile_residuals[key] == []:
            tile_mean[key] = np.nan
            unconnected_tile_ids.append(key)
        else:
            tile_mean[key] = np.mean(tile_residuals[key])

    return tile_mean, unconnected_tile_ids



def compute_residuals1(tilespecs_p, tilespecs_q, tile_residuals_p, tile_residuals_q, match):
    pts_p = np.array(match['matches']['p'])
    pts_q = np.array(match['matches']['q'])

    ts_p = next(ts for ts in tilespecs_p if ts.tileId == match['pId'])
    ts_q = next(ts for ts in tilespecs_q if ts.tileId == match['qId'])

    if ts_p is not None and ts_q is not None:
        t_p = ts_p.tforms[-1].tform(pts_p.T)
        t_q = ts_q.tforms[-1].tform(pts_q.T)

        all_pts = np.concatenate([t_p[1:,:], t_q[1:,:]], axis=1)
        res = np.sqrt(np.sum(np.power(all_pts[:,0:2]-all_pts[:,2:4],2),axis=1))

        tile_residuals_p[match['pId']] = np.append(tile_residuals_p[match['pId']], res)
        tile_residuals_q[match['qId']] = np.append(tile_residuals_q[match['qId']], res)



def compute_residuals(render, stack, matchCollectionOwner, matchCollection, Z, min_points=10):
    z_p = Z[0]
    z_q = Z[1]

    # get section ids for z_p and z_q
    pgroupId = render.run(renderapi.stack.get_sectionId_for_z, stack, z_p)
    qgroupId = render.run(renderapi.stack.get_sectionId_for_z, stack, z_q)

    print z_p, pgroupId

    allmatches = render.run(
                    renderapi.pointmatch.get_matches_from_group_to_group,
                    matchCollection,
                    pgroupId,
                    qgroupId,
                    owner=matchCollectionOwner)

    print allmatches[0]
    # tile specs are needed to extract the transformations
    tilespecs_p = render.run(renderapi.tilespec.get_tile_specs_from_z, stack, z_p)

    # get tile ids for each z
    # avoid reading the same tilespecs if z_p == z_q
    tilespecs_q = []
    tile_ids_q = []
    tile_ids_p = get_tile_ids(tilespecs_p)
    if z_p == z_q:
        tile_ids_q = tile_ids_p
        tilespecs_q = tilespecs_p
    else:
        tilespecs_q = render.run(renderapi.tilespec.get_tile_specs_from_z, stack, z_q)
        tile_ids_q = get_tile_ids(tilespecs_q)

    # filter and compute transformed points
    transformed_pts_p = np.zeros((1,2))
    transformed_pts_q = np.zeros((1,2))
    #transformed_pts_p = []
    #transformed_pts_q = []

    # initialize tile based residuals
    tile_residuals_p = {key: np.empty((0,1)) for key in tile_ids_p}
    tile_residuals_q = {key: np.empty((0,1)) for key in tile_ids_q}

    print "Computing residuals for every point match"
    #mypartial = partial(compute_residuals1, tilespecs_p, tilespecs_q, tile_residuals_p, tile_residuals_q)

    #with mp.ProcessingPool(20) as pool:
    #    pool.map(mypartial, allmatches)

    for i, match in enumerate(allmatches):
      #print i
      pts_p = np.array(match['matches']['p'])
      pts_q = np.array(match['matches']['q'])

      if pts_p.shape[1] < min_points:
          continue

      '''
      # get the transformation list for this point match
      try:
          ts_p = next(ts for ts in tilespecs_p if ts.tileId == match['pId'])
          ts_q = next(ts for ts in tilespecs_q if ts.tileId == match['qId'])
          print "H ", ts_p, ts_q
      except:
          try:
              ts_p = next(ts for ts in tilespecs_p if ts.tileId == match['qId'])
              ts_q = next(ts for ts in tilespecs_q if ts.tileId == match['pId'])
          except:
              ts_p = None
              ts_q = None
      '''

      ts_p = next(ts for ts in tilespecs_p if ts.tileId == match['pId'])
      ts_q = next(ts for ts in tilespecs_q if ts.tileId == match['qId'])

      if ts_p is not None and ts_q is not None:
          t_p = ts_p.tforms[-1].tform(pts_p.T)
          t_q = ts_q.tforms[-1].tform(pts_q.T)

          #transformed_pts_p = np.append(transformed_pts_p, t_p, axis=0)
          #transformed_pts_q = np.append(transformed_pts_q, t_q, axis=0)

          # tile based residual
          all_pts = np.concatenate([t_p[1:,:], t_q[1:,:]], axis=1)

          res = np.sqrt(np.sum(np.power(all_pts[:,0:2]-all_pts[:,2:4],2),axis=1))

          tile_residuals_p[match['pId']] = np.append(tile_residuals_p[match['pId']], res)
          tile_residuals_q[match['qId']] = np.append(tile_residuals_q[match['qId']], res)

    #transformed_pts_p = transformed_pts_p[1:,:]
    #transformed_pts_q = transformed_pts_q[1:,:]
    #all_points = np.concatenate([transformed_pts_p, transformed_pts_q], axis=1)
    #dv=np.sqrt(np.sum(np.power(all_points[:,0:2]-all_points[:,2:4],2),axis=1))
    return tile_residuals_p, tile_residuals_q



def calculate_montage_point_match_residuals(render, stack, matchCollectionOwner, matchCollection, zstart=None, zend=None, cutoff=3, normalizer=1):

    zvalues = render.run(renderapi.stack.get_z_values_for_stack, stack)
    if zstart is None:
        zstart = min(zvalues)
    elif zstart < min(zvalues):
        zstart = min(zvalues)

    if zend is None:
        zend = max(zvalues)
    elif zend > max(zvalues):
        zend = max(zvalues)

    # match group ids with their corresponding z
    # for EM this is straight forward as section(group) ids are z values
    # for AT the section Ids are different
    #section_ids = [render.run(renderapi.stack.get_sectionId_for_z, stack, z) for z in range(zstart, zend+1)]

    # for every z compute the residuals and compute overall diagnostics information for each section
    stats = {}
    for z in range(zstart, zend+1):
        zz = z - zstart
        stats[zz] = {}
        stats[zz]['group_mean'] = []
        stats[zz]['group_median'] = []
        stats[zz]['group_variance'] = []
        stats[zz]['outlier_counts'] = []
        stats[zz]['outlier_percent'] = []
        stats[zz]['unconnected_count'] = []
        stats[zz]['tile_residual_mean'] = []
        stats[zz]['tile_ids'] = []
        stats[zz]['outlier_tile_ids'] = []
        stats[zz]['unconnected_tile_ids'] = []

    opts = {'normalizer':1, 'method': 'std', 'cutoff':cutoff, 'greater_than':True}

    mypartial = partial(compute_residuals, render, stack, matchCollectionOwner, matchCollection)
    zvalues = range(zstart, zend+1)

    #with mp.ProcessingPool(5) as pool:
    #    tile_residuals_p, tile_residuals_q = pool.map(mypartial, it.izip(zvalues, zvalues))

    tile_residuals_p = []
    tile_residuals_q = []
    for z in range(zstart, zend+1):
        trp, trq = compute_residuals(render, stack, matchCollectionOwner, matchCollection, [z,z])
        print len(trp)
        tile_residuals_p.append(trp)
        tile_residuals_q.append(trq)

        #print tile_residuals_p[z-zstart]
        stats[z-zstart]['tile_residual_mean'], stats[z-zstart]['unconnected_tile_ids'] = compute_mean_tile_residuals(tile_residuals_p[z-zstart])
        #print stats[z-zstart]['tile_residual_mean']
        stats[z-zstart]['z'] = z
        stats[z-zstart]['group_mean'] = np.nanmean(stats[z-zstart]['tile_residual_mean'].values())
        #print stats[z-zstart]['tile_residual_mean']
        #stats[z-zstart]['group_variance'] = np.nanvar(stats[z-zstart]['tile_residual_mean'].values(), ddof = normalizer)
        #stats[z-zstart]['group_median'] = np.nanmedian(stats[z-zstart]['tile_residual_mean'].values())
        #stats[z-zstart]['group_std_deviation'] = np.sqrt(stats[z-zstart]['group_variance'])

    return stats, zvalues





if __name__ == "__main__":
    render = renderapi.Render(**example['render'])

    t1 = time.time()
    out, zvalues = calculate_montage_point_match_residuals(
            render,
            example['input_stack'],
            example['matchCollectionOwner'],
            example['matchCollection'],
            zstart=example['zstart'],
            zend=example['zend'])

    t2 = time.time()
    print "Time taken for computation %s seconds"%(str(t2-t1))

    print len(out)
    transformed_positions = calculate_section_maps(render, example['input_stack'], example['zstart'])

    plot_section_maps(transformed_positions, out[0]['tile_residual_mean'].values())
    #np.savetxt('trans_positions.txt', transformed_positions)

    #group_mean = []
    #for z in zvalues:
    #    group_mean.append(out[z-example['zstart']]['group_mean'])
    #plot_stack_statistics(group_mean, zvalues)
