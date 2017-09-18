# -*- coding: utf-8 -*-
"""
Created on Fri Jul 21 10:10:00 2017

@author: benjaminp
"""
import bokeh.plotting as bkp
import diagnostics
import numpy as np

from bokeh.io import output_file, show
from bokeh.layouts import column
from bokeh.plotting import figure, ColumnDataSource
from bokeh.models import Title
from bokeh.models import HoverTool

plot_cs_pm = False

montaged = {
        "render": {
            "owner":"gayathri",
            "project":"EM_Phase1",
            "host":"http://10.128.124.14",
            "port":8999,
            "client_scripts":"/data/nc-em2/gayathrim/Janelia_Pipeline/render/render-ws-java-client/src/main/scripts"
        },
        
        "matchSource":'Phase1Data_2366_2415_montagePM',
        'sourceStack':'Phase1Data_2366_2415_Montage'
}

simple_mat = np.load('save_simple_mat.npy')
avg_simple_mat = np.mean(simple_mat, axis = 1)
save_all_medians = np.load('save_all_medians.npy')

output_file('test.html')

def simple_plot(title = 'title', x_label = 'Section (z)', y_label = 'value', width = 800, height = 400, size = 15, tools = []):
    
    if tools != []:
        p = figure(plot_height = height, plot_width = width, tools = tools)
    else:
        p = figure(plot_height = height, plot_width = width)
    p.title.text = title
    p.title.align = "center"
    p.title.text_font_size = "25px"
    p.add_layout(Title(text=x_label, align="center"), "below")
    p.add_layout(Title(text=y_label, align="center"), "left")

    return p

#==============================================================================
# Cross-section PM Plotting
#==============================================================================

if plot_cs_pm:
    simple_mat = save_simple_mat

    rng = range(2268, 3481+1)
    lvl_1_resids = simple_mat[:,0:2]
    lvl_2_resids = simple_mat[:,2:4]
    means_lvl_1 = np.mean(lvl_1_resids, axis = 1)
    means_lvl_2 = np.mean(lvl_2_resids, axis = 1)
    means = np.column_stack((means_lvl_1, means_lvl_2))
    
    means_means = np.mean()
    
    cs_pm_source = ColumnDataSource(data = dict(xs = rng,
                                                lvl_1_means = means_lvl_1,
                                                lvl_2_means = means_lvl_2,
                                                lvl_1_p = simple_mat[:,0],
                                                lvl_1_q = simple_mat[:,1],
                                                lvl_2_p = simple_mat[:,2],
                                                lvl_2_q = simple_mat[:,3],
                                                ))
    
    cs_hov = HoverTool(tooltips = [
            ('section', '@xs'),
            ('1 level resid', '@lvl_1_means'),
            ('2 level resid', '@lvl_2_means'),
            ])
    
    
    cs_pm_fig = simple_plot('Cross-section residuals', y_label = 'Residual', tools = [cs_hov])
    
    cs_pm_fig.circle('xs', 'lvl_1_means', source = cs_pm_source, color = 'red', legend = '1 section step', )
    cs_pm_fig.circle('xs', 'lvl_2_means', source = cs_pm_source, legend = '2 section step', )
    cs_pm_fig.legend.click_policy = 'hide'
    #cs_pm_fig.circle('xs', 'lvl_1_p', color = 'red', legend = '1 section step p', source = cs_pm_source)
    #cs_pm_fig.circle('xs', 'lvl_2_p', color = 'yellow', legend = '2 section step p', source = cs_pm_source)
    #cs_pm_fig.circle('xs', 'lvl_1_q', color = 'blue', legend = '1 section step q', source = cs_pm_source)
    #cs_pm_fig.circle('xs', 'lvl_2_q', color = 'green', legend = '2 section step q', source = cs_pm_source)
    
    show(cs_pm_fig)

#==============================================================================
# Montage PM plotting
#==============================================================================
starts = range(2266,3416,50)
ends = range(2315,3465,50)
# variance, mean, median, outlier counts, outlier percents, unconnected counts
start = 2266
end = 3415

#start = 2268
#end = 3481
all_sections = range(start,end + 1)


current_medians = save_all_medians
#compare_simple_mat = avg_simple_mat[:1148]

mean_fig = simple_plot('Mean montage PM residuals for section %i - %i' %(start, end), y_label = 'Residual')
mean_fig.circle(all_sections, current_medians)
show(mean_fig)
#
#compare_cs_mont = simple_plot('Compare', x_label = 'montage', y_label = 'cross-section')
#compare_cs_mont.circle(current_medians, compare_simple_mat)
#show(compare_cs_mont)
#corr = np.corrcoef(compare_simple_mat[0:300], current_medians[0:300])