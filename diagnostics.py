# -*- coding: utf-8 -*-
"""
Created on Fri Jul 07 09:09:30 2017

7/24/2017 update


@author: Benjamin Pedigo
"""


import numpy as np
import renderapi as ren
import warnings
import bokeh.plotting as bkp
warnings.simplefilter(action='ignore', category=FutureWarning)

def get_tile_ids(tiles_list):
    """Return an array of unicode tile ids given a list of tilespec objects."""
    id_array = np.empty(len(tiles_list),dtype = '<U64')
    for i, name in enumerate(tiles_list):
        #print name.tileId
        id_array[i] = name.tileId
    return id_array

# return array of unicode tile ids given list of point matches
def get_tile_ids_from_matches(matches_list, key, unique = True):
    """Return an array of unicode tile ids given a list of point match dictionaries.

    matches_list -- the list of point match dicts
    key -- what key to find in the PM dicts, usually 'pId' or 'qId'
    unique -- bool; whether or not to remove duplicates from the output array
    """
    id_array = np.empty(len(matches_list),dtype = '<U64')
    for i, match in enumerate(matches_list):
        id_array[i] = match[key]
    if unique:
        id_array = np.unique(id_array)
    return id_array


def tform_poly2d_asAffine(points, pr):
    '''Returns array of x, y points after applying an transformation specified
    by pr (parameters of transform). The transformation is done only using the
    affine terms from a poly2d transform object.

    points -- array of x,y paired points to be transformed
    pr -- parameters of the transformation from 2d polynomial transform obj
    '''
    points_wZ = np.column_stack((points, np.ones(len(points[:,0]))))

    affine_approx_Tmatrix = np.array([[pr[0,1], pr[1,1], 0], [pr[0,2], pr[1,2], 0], [pr[0,0], pr[1,0], 1]])

    affine_shifted_points = np.dot(points_wZ, affine_approx_Tmatrix)

    return affine_shifted_points

def stats_and_outliers(data, ids, cutoff = 1, method = 'std', greater_than = True, normalizer = 1, **kwargs):
    '''Returns dict of various statistics and information about outliers given
    an input dataset. The input data in this case is usually a list of tile
    ratios (area and perimeter, for example) or list of tile residuals (from
    point matches)

    cutoff -- the value to compare data to to determine whether or not it is
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

def polyarea(x,y):
    '''Returns the area of a polygon (usually 4 corners here but not necessarily
    based on input verticies'''

    return 0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))

def perimeter(x,y):
    '''Returns the perimeter of a polygon when given ordered lists of x,y verts.
    Must be in an order that does not cross through the center of a shape that
    the user is specifying.
    '''
    per = 0
    j = 0 + 1j
    for i in range(len(x)-1):
        per += abs(x[i+1] - x[i] + (y[i+1] - y[i])*j)

    per += abs(x[0]-x[3] + (y[0] - y[3])*j)
    return per

def plot_by_group_id(data, group_ids):
    bkp.output_file('test.html')
    p = bkp.figure()
    p.circle(group_ids, data, size = 20)
    bkp.show(p)

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


def calculate_montage_pm_residuals(mont_pm, z_start = None, z_end = None, normalizer = 1, method = 'std', cutoff = 3, greater_than = True, verbose = True, **kwargs):
    '''Returns summary dict of information about the point match residual
    statistics for a specified Render point match collection.

    z_start/z_end -- z values to calculate statistics on
    normalizer/method/cutoff/greater_that -- see function stats_and_outliers
    verbose -- whether or not to output progress through Z sections

    Output dictionary
    'means' -- array with mean pm residual for each section. Calculated as the
    mean of each tile's mean residual
    'medians' -- array with median pm residual for each section. Median of all
    tile's mean residual
    'means for each tile' -- list of arrays, listed by z section, each array
    contains the mean residual for each tile
    'outlier counts' -- array of number of outliers detected for each section
    'outlier percent' -- above in percentage
    'outlier tile ids' -- array of unicode tile ids that were flagged
    'outlier tile idx' -- index of those tiles in the tile ids list
    'tile ids' -- list of arrays, arrays contain ids for each tile in the sect.
    'unconnected tile count' -- number of tiles with no point matches within
    the tile
    'unconnected tile ids' -- array of unicode tile ids for the above
    'unconnected tile idx' -- indeces of the above in tile ids
    'values' -- list of lists of lists... the first dimensions is a list by
    section, the next is by tile, and the last is a list containing the
    residuals for each point match for that tile
    'variances' -- array of variances of residuals for each section
    '''

    r = ren.Render(**mont_pm['render'])
    group_ids = ren.pointmatch.get_match_groupIds_from_only(mont_pm['matchSource'], owner=mont_pm['matchCollectionOwner'], render=r)
    group_ids = [float(i) for i in group_ids]

    if method == 'fixed':
            cutoff += -1 # because code assumes ratio input
    # repack opts
    opts = {'normalizer':normalizer, 'method': method, 'cutoff':cutoff, 'greater_than':greater_than}


    if z_start == None or z_end == None:
        z_start = min(group_ids)
        z_end = max(group_ids)


    z_first_idx = [idx for idx,ids in enumerate(group_ids) if ids == float(z_start)][0]
    z_len = int(z_end - z_start) + 1


    groups_mean = np.zeros(z_len)
    groups_median = np.zeros(z_len)
    groups_variance = np.zeros(z_len)
    groups_outlier_counts = np.zeros(z_len)
    groups_outlier_percent = np.zeros(z_len)
    groups_unconnected_count = np.zeros(z_len)

    # pre allocating lists for parallelization?
    groups_all_tiles_ids = [None] * z_len
    groups_outlier_tiles_ids = [None] * z_len
    groups_outlier_tiles_idx = [None] * z_len
    groups_all_tiles_means = [None] * z_len
    groups_unconnected_tiles_ids = [None] * z_len
    groups_unconnected_tiles_idx = [None] * z_len
    groups_values = [None] * z_len

    for j in range(z_len): #
        group_id = group_ids[j + z_first_idx]

        if verbose: print 'Calculating montage residuals for Z = %i' %(int(float(group_id)))


        matches = ren.pointmatch.get_matches_within_group(mont_pm['matchSource'],group_id, owner=mont_pm['matchCollectionOwner'], render=r)

        tiles = ren.tilespec.get_tile_specs_from_z(mont_pm['sourceStack'], float(group_id), render = r)

        tiles_id = get_tile_ids(tiles)

        tile_residuals = [[] for x in range(len(tiles_id))]

        for x in range(len(matches)): #len(matches)

            match = matches[x]

            residual, p_idx, q_idx = calculate_residual_for_match(match,tiles_id,tiles)
            if residual =='flag':
                continue

            tile_residuals[p_idx].append(residual)
            tile_residuals[q_idx].append(residual)


        tile_mean = np.zeros(len(tile_residuals))
        unconnected_tile_idx = []
        for i in range(len(tile_residuals)):
            if tile_residuals[i] == []:
                tile_mean[i] = np.nan
                unconnected_tile_idx.append(i)
            else:
                tile_mean[i] = np.mean(tile_residuals[i]) # mean for a single tile

        unconnected_count = len(unconnected_tile_idx)

        unconnected_tile_ids = tiles_id[unconnected_tile_idx]

        group_summary = stats_and_outliers(tile_mean, tiles_id, **opts)

        #construct output
        groups_mean[j] = group_summary['mean']
        groups_median[j] = group_summary['median']
        groups_variance[j] = group_summary['variance']
        groups_outlier_counts[j] = group_summary['outlier count']
        groups_outlier_percent[j] = group_summary['outlier percent']
        groups_unconnected_count[j] = unconnected_count
        groups_all_tiles_ids[j] = tiles_id
        groups_all_tiles_means[j] = tile_mean
        groups_outlier_tiles_ids[j] = group_summary['outlier ids']
        groups_outlier_tiles_idx[j] = group_summary['outlier idx']
        groups_unconnected_tiles_ids[j] = unconnected_tile_ids
        groups_unconnected_tiles_idx[j] = unconnected_tile_idx
        groups_values[j] = tile_residuals

    output = {'means':groups_mean, 'medians':groups_median,
              'variances':groups_variance, 'outlier counts':groups_outlier_counts,
              'outlier percent':groups_outlier_percent, 'unconnected tile count':groups_unconnected_count,
              'tile ids': groups_all_tiles_ids, 'means for each tile':groups_all_tiles_means,
              'outlier tile ids':groups_outlier_tiles_ids, 'outlier tile idx':groups_outlier_tiles_idx,
              'unconnected tile ids':groups_unconnected_tiles_ids, 'unconnected tile idx':groups_unconnected_tiles_idx, 'values':groups_values}

    #plot_by_group_id(groups_mean, group_ids)

    return output


def calculate_area_perimeter_ratios(raw, montaged, z_start = None, z_end = None,
                                    normalizer = 1, method = 'std', cutoff = 3,
                                    greater_than = True, verbose = True, *args, **kwargs):
    '''Returns a dictionary of various information and stats on the area and
    perimeter ratios (ratio of deformed (target) tile over undeformed)

    Note: uses linear approximation of tile deformation even if polynomial
    transform was found. (same as EM_aligner)

    Note: the 'source' tile, or denomenator in the ratio, is assumed to be
    rectangular. (same as EM_aligner)

    normalizer/method/cutoff/greater_than -- see function stats_and_outliers
    verbose -- whether to output progress by z section


    Output dict
    'area' -- another dictionary with area stats/info
    'perimeter' -- another dictionary with perimeter stats/info

    Within each of the above dictionaries:
        'means' -- array containing mean ratio for each section
        'medians' -- array containing median ration for each section
        'outlier count' -- number of outliers specified by chosen method, array
        with a count for each section
        'outlier ids' -- indices of above in tile
        'outlier percent' -- the above count as a percentage of tiles in section
        'ratios' -- list of arrays (one array per section) of ratios for each
        tile
        'target values' -- list of arrays (one array per secion) where elements
        of arrays are non-ratio area or perimeter values for each tile
        'variances' -- array of variances (of the ratios) for each section
    '''

    # repack opts - couldn't find a better way to do this but there may be one
    opts = {'normalizer':normalizer, 'method': method, 'cutoff':cutoff, 'greater_than':greater_than}

    r_raw = ren.Render(**raw['render'])
    r_mont = ren.Render(**montaged['render'])

    if z_start == None or z_end == None:
        rawBounds = ren.stack.get_stack_bounds(raw['sourceStack'], render=r_raw)
        montageBounds = ren.stack.get_stack_bounds(montaged['sourceStack'], render=r_mont)

        z_start = int(max(rawBounds['minZ'], montageBounds['minZ']))
        z_end = int(min(rawBounds['maxZ'], montageBounds['maxZ']))

    z_list = range(z_start, z_end + 1)
    Zs = len(z_list)

    # lists
    Zs_area_target = [None] * Zs
    Zs_area_ratio = [None] * Zs
    Zs_area_outlier_ids = [None] * Zs
    Zs_area_outlier_idx = [None] * Zs

    Zs_perimeter_target = [None] * Zs
    Zs_perimeter_ratio = [None] * Zs
    Zs_perimeter_outlier_ids = [None] * Zs
    Zs_perimeter_outlier_idx = [None] * Zs


    # arrays
    Zs_area_mean = np.zeros(Zs)
    Zs_area_median = np.zeros(Zs)
    Zs_area_variance = np.zeros(Zs)
    Zs_area_outlier_count = np.zeros(Zs)
    Zs_area_outlier_percent = np.zeros(Zs)

    Zs_perimeter_mean = np.zeros(Zs)
    Zs_perimeter_median = np.zeros(Zs)
    Zs_perimeter_variance = np.zeros(Zs)
    Zs_perimeter_outlier_count = np.zeros(Zs)
    Zs_perimeter_outlier_percent = np.zeros(Zs)

    for z_index in range(Zs):

        if verbose:
            print('Calculating area and perimeter ratios for Z = %i' %(z_list[z_index]))

        tiles_source = ren.tilespec.get_tile_specs_from_z(raw['sourceStack'], z_list[z_index], render = r_raw)
        tiles_target = ren.tilespec.get_tile_specs_from_z(montaged['sourceStack'], z_list[z_index], render = r_mont)

        ids_source = get_tile_ids(tiles_source)
        ids_target = get_tile_ids(tiles_target)

        ids_source_sorted, ids_source_idx = np.unique(ids_source, return_index=True)
        is_source_in_target = np.in1d(ids_source, ids_target, assume_unique=True)
        ids_source_map = ids_source_idx[is_source_in_target]

        area_ratio = np.zeros(len(ids_target))
        area_target = np.zeros(len(ids_target))
        perimeter_target = np.zeros(len(ids_target))
        perimeter_ratio = np.zeros(len(ids_target))

        # loop over all tiles in a slice
        for tile_index in range(len(ids_target)):
            tile_target = tiles_target[tile_index]
            tile_source = tiles_source[ids_source_map[tile_index]]

            # assuming untransformed original tile
            corners = np.array([[0, 0], [tiles_target[0].width, 0], [tile_target.width,
                                tile_target.height],[0, tile_target.height]])

            # approximating as affine
            # note: assumes that the original tile was untransformed so affine_corners_source is not yet used
            if tile_target.tforms[1].className == 'mpicbg.trakem2.transform.PolynomialTransform2D': #and tile_source.tforms[1].className == 'mpicbg.trakem2.transform.PolynomialTransform2D':
                affine_corners = tform_poly2d_asAffine(corners, tile_target.tforms[1].params)
                #affine_corners_source = tform_poly2d_asAffine(corners, tile_source.tforms[1].params)
            elif tile_target.tforms[1].className == 'mpicbg.trakem2.transform.AffineModel2D': #and tile_source.tforms[1].className == 'mpicbg.trakem2.transform.AffineModel2D':
                affine_corners = tile_target.tforms[1].tform(corners)
                #affine_corners_source = tile_source.tforms[1].tform(corners)
            else:
                print 'Could not compatible transform 2D in second transform slot, no transformation taken'


            area_target[tile_index] = polyarea(affine_corners[:,0],affine_corners[:,1])
            area_ratio[tile_index] = area_target[tile_index] / (
                    tile_source.width * tile_source.height)


            perimeter_target[tile_index] = perimeter(affine_corners[:,0],affine_corners[:,1])
            perimeter_ratio[tile_index] = perimeter_target[tile_index]/perimeter(corners[:,0],corners[:,1])



        # load into data structures across Z
        Zs_area_target[z_index]=area_target
        Zs_area_ratio[z_index]=area_ratio


        area_stats_outliers = stats_and_outliers(area_ratio, ids_target, **opts)
        Zs_area_mean[z_index] = area_stats_outliers['mean']
        Zs_area_variance[z_index] = area_stats_outliers['variance']
        Zs_area_median[z_index] = area_stats_outliers['median']
        Zs_area_outlier_idx[z_index] = area_stats_outliers['outlier idx']
        Zs_area_outlier_ids[z_index] = area_stats_outliers['outlier ids']
        Zs_area_outlier_count[z_index] = area_stats_outliers['outlier count']
        Zs_area_outlier_percent[z_index] = area_stats_outliers['outlier percent']


        Zs_perimeter_target[z_index] = perimeter_target
        Zs_perimeter_ratio[z_index] = perimeter_ratio

        perimeter_stats_outliers = stats_and_outliers(perimeter_ratio, ids_target, **opts)
        Zs_perimeter_mean[z_index] = perimeter_stats_outliers['mean']
        Zs_perimeter_variance[z_index] = perimeter_stats_outliers['variance']
        Zs_perimeter_median[z_index] = perimeter_stats_outliers['median']
        Zs_perimeter_outlier_idx[z_index] = perimeter_stats_outliers['outlier idx']
        Zs_perimeter_outlier_ids[z_index] = perimeter_stats_outliers['outlier ids']
        Zs_perimeter_outlier_count[z_index] = perimeter_stats_outliers['outlier count']
        Zs_perimeter_outlier_percent[z_index] = perimeter_stats_outliers['outlier percent']


    # make output dictionary
    area_dict = {'target values':Zs_area_target, 'ratios':Zs_area_ratio, 'means':Zs_area_mean, 'variances':Zs_area_variance, 'medians':Zs_area_median, 'outlier count':Zs_area_outlier_count, 'outlier ids': Zs_area_outlier_ids, 'outlier percent':Zs_area_outlier_percent, 'outlier indecies':Zs_area_outlier_idx}
    perimeter_dict = {'target values':Zs_perimeter_target, 'ratios':Zs_perimeter_ratio, 'means':Zs_perimeter_mean, 'variances':Zs_perimeter_variance, 'medians':Zs_perimeter_median, 'outlier count':Zs_perimeter_outlier_count, 'outlier ids': Zs_perimeter_outlier_ids, 'outlier percent':Zs_perimeter_outlier_percent, 'outlier indecies':Zs_perimeter_outlier_idx}
    output = {'area':area_dict, 'perimeter':perimeter_dict}

    return output


def calculate_cross_sec_pm_residuals(fine_pm, z_start = None, z_end = None, num_cross_sections = 2, verbose = True, *args, **kwargs):
    '''Returns matrix (2d array) of residuals between sections.

    num_cross_sections -- default to 2; this specifies how many sections above
    to search for cross-section PMs.
    verbose -- bool; whether or not to output progress for each Z section


    Note: Currently cannot find cross-section point matches for the last 2 Z
    values in a stack.

    resid_mat --
    Residuals are calculated between a section (p) and the one above it (q).
    Each tile in the p section will have several point matches and therefore
    several residuals. The mean residual is taken for each tile in p and q
    (p/q_mean_residuals). Even though the residual for each point match is the
    same for the p or the q tile, the MEAN is taken by TILE, so the matrix
    below will not necessarily be symmetric. The MEDIAN is then taken of each
    of these lists, and stored in the appropriate slot in the matrix. In the
    matrix below, the section that was used for the mean/median is denoted with
    a capitol P/Q. For example, row 1 column 2: the point matches were between
    Z = 0 and Z = 1, but the means were taken for each tile in Z = 0.

    -------------------------------------------------------------------------
    |       nan       | P = 0 --> q = 1 | P = 0 --> q = 2 |       nan       |
    -------------------------------------------------------------------------
    | p = 0 --> Q = 1 |       nan       | P = 1 --> q = 2 | P = 1 --> q = 3 |
    -------------------------------------------------------------------------
    | p = 0 --> Q = 2 | p = 1 --> Q = 2 |       nan       | P = 2 --> q = 3 |
    -------------------------------------------------------------------------
    |       nan       | p = 1 --> Q = 3 | p = 2 --> Q = 3 |       nan       |
    -------------------------------------------------------------------------

    simple_resid_mat -- a simpler version of the above
    each row is a z-section with 4 columns:
    -------------------------------------------------------------------------
    | P = 0 --> q = 1 | p = 0 --> Q = 1 | P = 0 --> q = 2 | p = 0 --> Q = 2 |
    -------------------------------------------------------------------------
    | P = 1 --> q = 1 | p = 1 --> Q = 1 | P = 1 --> q = 2 | p = 1 --> Q = 2 |
    -------------------------------------------------------------------------

    '''
    np.seterr(divide='ignore', invalid='ignore')

    r = ren.Render(**fine_pm['render'])

    group_ids = ren.pointmatch.get_match_groupIds_from_only(fine_pm['matchSource'], render=r)
    group_ids = [float(i) for i in group_ids]

    # currently, this doesn't work because no point matches are listed outside of the stack in what I am accessing from Render
    #group_ids.append(group_ids[-1] + 1)
    #group_ids.append(group_ids[-1] + 1)

    if z_start == None or z_end == None:
        z_start = min(group_ids)
        z_end = max(group_ids) - num_cross_sections # exclusive here


    z_first_idx = [idx for idx,ids in enumerate(group_ids) if ids == float(z_start)][0]
    z_len = int(z_end - z_start) + 1

    resid_mat = np.empty([z_len, z_len])
    resid_mat[:] = np.NAN
    simple_resid_out = np.zeros([z_len, num_cross_sections*2])

    for z_idx in range(z_len):

        if verbose:
            print 'Calculating cross section residuals for Z = %i' %(int(float(group_ids[z_idx + z_first_idx])))

        p_tiles = ren.tilespec.get_tile_specs_from_z(fine_pm['sourceStack'], group_ids[z_idx +z_first_idx], render=r)
        p_tiles_id = get_tile_ids((p_tiles))

        for step_idx in range(num_cross_sections):
            q_tiles = ren.tilespec.get_tile_specs_from_z(fine_pm['sourceStack'], group_ids[z_idx + step_idx + 1 + z_first_idx], render=r)
            q_tiles_id = get_tile_ids(q_tiles)

            matches = ren.pointmatch.get_matches_from_group_to_group(fine_pm['matchSource'], group_ids[z_idx + z_first_idx], group_ids[z_idx + step_idx + 1 + z_first_idx], render = r)
            # can be used to check if tile ids found in the match collection correspond to those in the stack (usually no)
#            p_tiles_id_from_matches = get_tile_ids_from_matches(matches, 'pId')
#            q_tiles_id_from_matches = get_tile_ids_from_matches(matches, 'qId')


            p_residuals = np.zeros(len(p_tiles))
            q_residuals = np.zeros(len(q_tiles))

            p_counts = np.zeros(len(p_tiles))
            q_counts = np.zeros(len(q_tiles))

            for match in matches:

                residual, p_idx, q_idx = calculate_residual_for_match(match,p_tiles_id, p_tiles, q_tiles_id, q_tiles)
                if residual == 'flag':
#                    print 'Could not find tile id: %s or %s in corresponding tile spec' %(str(match['pId']), str(match['qId']))
#                    print 'P is %i, Q is %i' %(z_idx, step_idx)
#                    print ''
                    continue

                p_residuals[p_idx] += residual
                q_residuals[q_idx] += residual

                p_counts[p_idx] += 1
                q_counts[q_idx] += 1


            # Note: using nanmedian becuase some entries in p_counts end up with 0,
            # unsure why this does not occur for Matlab code but this seems to
            # end up with the same result
            p_mean_residuals = np.divide(p_residuals, p_counts)
            p_median_total_residuals = np.nanmedian(p_mean_residuals)


            q_mean_residuals = np.divide(q_residuals, q_counts)
            q_median_total_residuals = np.nanmedian(q_mean_residuals)

            if z_idx+step_idx < z_len - 1:
                resid_mat[z_idx, z_idx + step_idx + 1] = p_median_total_residuals
                resid_mat[z_idx + step_idx + 1, z_idx] = q_median_total_residuals

            simple_resid_out[z_idx, 0 + step_idx*2] = p_median_total_residuals
            simple_resid_out[z_idx, 1 + step_idx*2] = q_median_total_residuals
    return resid_mat, simple_resid_out
