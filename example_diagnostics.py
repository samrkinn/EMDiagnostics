# -*- coding: utf-8 -*-
"""
Created on Wed Jul 12 15:48:31 2017

@author: Benjamin Pedigo
"""

import diagnostics
import time
import json

raw = {
    "render": {
        "owner":"gayathri",
        "project":"EM_Phase1",
        "host":"http://10.128.124.14",
        "port":8999,
        "client_scripts":"/data/nc-em2/gayathrim/Janelia_Pipeline/render/render-ws-java-client/src/main/scripts"
    },

    "sourceStack":"Phase1Data_2316_2365",
}

montaged = {
    "render": {
        "owner":"gayathri",
        "project":"MM2",
        "host":"http://10.128.124.14",
        "port":8998,
        "client_scripts":"/data/nc-em2/gayathrim/Janelia_Pipeline/render/render-ws-java-client/src/main/scripts"
    },
    "matchCollectionOwner":"gayathri_MM2",
    "matchSource":"mm2_acquire_8bit_montage",
    'sourceStack':"mm2_acquire_8bit_Montage"
}

fine_pm = {
    "render": {
        "owner":"gayathri",
        "project":"EM_Phase1_Fine",
        "host":"http://10.128.124.14",
        "port":8999,
        "client_scripts":"/data/nc-em2/gayathrim/Janelia_Pipeline/render/render-ws-java-client/src/main/scripts"
    },

    "matchSource":"FineAlign_2316_2365_AFF_CONS_Scale_060",
    'sourceStack':"FineAlign_fusion1"
}

time1 = time.time()
area_perimeter_opts = {'cutoff': 0.01, 'method': 'fixed', 'greater_than': False, 'normalizer':1, 'verbose':True}
A_P_stats = diagnostics.calculate_area_perimeter_ratios(raw, montaged, **area_perimeter_opts)
print 'Area and perimeter ratios took ' + str(time.time() - time1) + ' seconds'
print ''

time2 = time.time()
montage_pm_opts = {'greater_than':False, 'method':'std', 'cutoff':3, 'normalizer':1, 'verbose':True}
montage_residuals = diagnostics.calculate_montage_pm_residuals(montaged, z_start=1015, z_end=1118, **montage_pm_opts)
print 'Montage point match residuals took  ' + str(time.time() - time2)  + ' seconds'
print ''

with open('montage_mm2_stats.json', 'w') as f:
    json.dump(list(montage_residuals['means']), f, indent=5)

#time3 = time.time()
#cross_section_pm_opts = {}
#cross_section_residuals = diagnostics.calculate_cross_sec_pm_residuals(fine_pm)
#print 'Cross-section point match residuals took  ' + str(time.time() - time3)  + ' seconds'
#print ''

print 'Total analysis took ' + str(time.time() - time1) + ' seconds'
