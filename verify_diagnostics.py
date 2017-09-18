# -*- coding: utf-8 -*-
"""
Created on Wed Jul 12 15:48:31 2017

@author: Benjamin Pedigo
"""

import diagnostics
import time
import scipy.io as sio
import numpy as np
import bokeh as bk

starts = range(2266,3416,50)
ends = range(2315,3465,50)

raw = {
    "render": {
        "owner":"gayathri",
        "project":"EM_Phase1",
        "host":"http://10.128.124.14",
        "port":8999,
        "client_scripts":"/data/nc-em2/gayathrim/Janelia_Pipeline/render/render-ws-java-client/src/main/scripts"
    },
    
    "sourceStack":"Phase1Data_2266_2315",
}
    
montaged = {
    "render": {
        "owner":"gayathri",
        "project":"EM_Phase1",
        "host":"http://10.128.124.14",
        "port":8999,
        "client_scripts":"/data/nc-em2/gayathrim/Janelia_Pipeline/render/render-ws-java-client/src/main/scripts"
    },
    
    "matchSource":'Phase1Data_2266_2315_montagePM',
    'sourceStack':'Phase1Data_2266_2315_Montage'
}

fine_pm = {
    "render": {
        "owner":"gayathri",
        "project":"EM_Phase1_Fine",
        "host":"http://10.128.124.14",
        "port":8999,
        "client_scripts":"/data/nc-em2/gayathrim/Janelia_Pipeline/render/render-ws-java-client/src/main/scripts"
    },
    
    "matchSource":"FineAlign_fullStack",
    'sourceStack':"FineAlign_fusion1"
}

compare_to_matlab = False
test_A_P = False
test_mont_pm = True
test_cs_pm = False

tol = 1e-12

raw = {
    "render": {
        "owner":"gayathri",
        "project":"EM_Phase1",
        "host":"http://10.128.124.14",
        "port":8999,
        "client_scripts":"/data/nc-em2/gayathrim/Janelia_Pipeline/render/render-ws-java-client/src/main/scripts"
    },
    
    "sourceStack":'',
}
    
montaged = {
    "render": {
        "owner":"gayathri",
        "project":"EM_Phase1",
        "host":"http://10.128.124.14",
        "port":8999,
        "client_scripts":"/data/nc-em2/gayathrim/Janelia_Pipeline/render/render-ws-java-client/src/main/scripts"
    },
    
    'sourceStack':''
}

verbose = True
total_time = time.time()
if test_A_P:
    AP_output_by_stack = [None] * len(ends)
    check_area_whole_stack = [None] * len(ends)
    check_perimeter_whole_stack = [None] * len(ends)
    area_perimeter_opts = {'cutoff': 0.04, 'method': 'fixed', 'greater_than': False, 'normalizer':1, 'verbose':verbose}
    
    for i in range(len(ends)):
        if verbose: print 'Calculating area and perimeter ratios for stack %i - %i' %(starts[i], ends[i])
        raw['sourceStack'] = "Phase1Data_%i_%i"%(starts[i], ends[i])
        montaged['sourceStack'] = 'Phase1Data_%i_%i_Montage'%(starts[i], ends[i])
        
        time1 = time.time()
        A_P_stats = diagnostics.calculate_area_perimeter_ratios(raw, montaged, starts[i],ends[i], **area_perimeter_opts)
        AP_output_by_stack[i] = A_P_stats
        if verbose: print 'Area and perimeter ratios took ' + str(time.time() - time1) + ' seconds'
        
        if compare_to_matlab:
            if verbose: print 'Equal to Matlab output for stack %i - %i:' %(starts[i], ends[i])
            
            a_file = 'a_output_%i_%i.mat' %(starts[i],ends[i])
            a_matlab = sio.loadmat(a_file)
            check_area_means = np.allclose(a_matlab['mean'].T, A_P_stats['area']['means'], atol= tol, equal_nan = True)
            check_area_medians = np.allclose(a_matlab['median'].T, A_P_stats['area']['medians'], atol= tol, equal_nan = True)
            check_area_variances = np.allclose(a_matlab['variance'].T, A_P_stats['area']['variances'], atol= tol, equal_nan = True)
            check_area_outlier_count = np.allclose(a_matlab['outlier_count'].T, A_P_stats['area']['outlier count'], atol= tol, equal_nan = True)
            check_area_all = [check_area_means, check_area_medians, check_area_variances, check_area_outlier_count]
            
            p_file = 'p_output_%i_%i.mat' %(starts[i],ends[i])
            p_matlab = sio.loadmat(p_file)
            check_perimeter_means = np.allclose(p_matlab['mean'].T, A_P_stats['perimeter']['means'], atol= tol, equal_nan = True)
            check_perimeter_medians = np.allclose(p_matlab['median'].T, A_P_stats['perimeter']['medians'], atol= tol, equal_nan = True)
            check_perimeter_variances = np.allclose(p_matlab['variance'].T, A_P_stats['perimeter']['variances'], atol= tol, equal_nan = True)
            check_perimeter_outlier_count = np.allclose(p_matlab['outlier_count'].T, A_P_stats['perimeter']['outlier count'], atol= tol, equal_nan = True)
            check_perimeter_all = [check_perimeter_means, check_perimeter_medians, check_perimeter_variances, check_perimeter_outlier_count]
            
            if np.all(check_area_all) and np.all(check_perimeter_all):
                check_area_whole_stack[i] = True 
                check_perimeter_whole_stack[i] = True 
                if verbose: print True
            else:
                check_area_whole_stack[i] = False
                check_perimeter_whole_stack[i] = False
                if verbose: print False
                
                
                
###############################################################################


tiles = sio.loadmat('tile_id')['tile_id']
tiles = np.ravel(tiles)
mat_tiles_id = np.hstack(tiles)

all_medians = []
if test_mont_pm:
    
    check_mont_pm_whole_stack = [None] * len(ends)

    montaged = {
        "render": {
            "owner":"gayathri",
            "project":"EM_Phase1",
            "host":"http://10.128.124.14",
            "port":8999,
            "client_scripts":"/data/nc-em2/gayathrim/Janelia_Pipeline/render/render-ws-java-client/src/main/scripts"
        },
        
        "matchSource":'Phase1Data_2266_2315_montagePM',
        'sourceStack':'Phase1Data_2266_2315_Montage'
    }
        
    montage_pm_opts = {'greater_than':True, 'method':'fixed', 'cutoff':30, 'normalizer':1, 'verbose':False}

    for i in range(len(ends)):
        montaged['sourceStack'] = 'Phase1Data_%i_%i_Montage'%(starts[i], ends[i])
        montaged['matchSource'] = 'Phase1Data_%i_%i_montagePM' %(starts[i], ends[i])
        
        
        time2 = time.time()
        if verbose: print 'Calculating montage point match residuals for stack %i - %i' %(starts[i], ends[i])
    
        montage_residuals = diagnostics.calculate_montage_pm_residuals(montaged, **montage_pm_opts)
        all_medians.append(montage_residuals['medians'])
        
        if verbose: print 'Montage point match residuals took  ' + str(time.time() - time2)  + ' seconds'
        
        
        if compare_to_matlab:
            pm_file = 'mont_pm_resid_%i_%i.mat' %(starts[i],ends[i])
            mont_pm_resid_matlab = sio.loadmat(pm_file)
            check_means = np.allclose(mont_pm_resid_matlab['mean_of_means'].T, montage_residuals['means'], atol= tol, equal_nan = True)
            check_medians = np.allclose(mont_pm_resid_matlab['median_of_means'].T, montage_residuals['medians'], atol = tol, equal_nan = True)
            check_variances = np.allclose(mont_pm_resid_matlab['variance_of_means'].T, montage_residuals['variances'], atol = tol, equal_nan = True)
            check_outlier_count = np.allclose(mont_pm_resid_matlab['outlier_count'].T, montage_residuals['outlier counts'], atol = tol, equal_nan = True)
            check_unconnected_count = np.allclose(mont_pm_resid_matlab['unconnected_count'].T, montage_residuals['unconnected tile count'], atol = tol, equal_nan = True)
            check_all = [check_means, check_medians, check_variances, check_outlier_count, check_unconnected_count]
            # note: not checking indecies at the moment for outliers or for unconnected tiles. also not comparing tile ids, or individual values for each tile
    
            if np.all(check_all):
                check_mont_pm_whole_stack[i] = True 
                if verbose: print True
            else:
                check_mont_pm_whole_stack[i] = False
                if verbose: print False
            

###############################################################################
if test_cs_pm:
    time3 = time.time()
    cross_section_pm_opts = {'verbose':verbose}
    cross_section_residuals, simple_mat = diagnostics.calculate_cross_sec_pm_residuals(fine_pm, 2268, 3481, **cross_section_pm_opts) # #2400, 2405 #2268 #3481 
    print 'Cross-section point match residuals took  ' + str(time.time() - time3)  + ' seconds'
    print ''
    
    
    if compare_to_matlab:
        matlab_cs_pm = sio.loadmat('cs_pm_matlab.mat')['resmat']
        check = np.isclose(matlab_cs_pm, cross_section_residuals, atol = tol, equal_nan = True)
        if np.all(check):
            print 'Matlab and Python cross-section residuals agree to at least %f precision' %(tol)
            check_cs_pm = True
        else:
            print 'Some element(s) of Matlab and Python cross-section residuals do not agree'
            no_match_idx = np.where(check == False)
            number_no_match = len(no_match_idx[0])
            percent_no_match = number_no_match/float(np.count_nonzero(~np.isnan(cross_section_residuals))) * 100
            check_cs_pm = False


if verbose: print 'Entire comparison took ' + str(time.time() - total_time)  + ' seconds'


