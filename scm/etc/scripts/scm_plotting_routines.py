#!/usr/bin/env python

# -*- coding: utf-8 -*-
"""
Created on Mon May 16 16:58:28 2016

@author: grantf
"""
from __future__ import print_function

import numpy as np
import matplotlib as mpl
mpl.use('PS')
import matplotlib.pyplot as plt

mpl.rcParams['savefig.dpi'] = 150
mpl.rcParams['legend.handlelength'] = 5
latex_labels = False

def calc_skill_score(data, experiment, R_0, k):
    sigma_data = np.std(data)
    sigma_experiment = np.std(experiment)
    R = np.corrcoef(x=data, y=experiment, bias=1)[0,1]

    return 4.0*(1.0 + R)**k/((sigma_data/sigma_experiment+sigma_experiment/sigma_data)**2*(1.0 + R_0)**k)
    #return np.log10(((sigma_data/sigma_experiment+sigma_experiment/sigma_data)**2*(1.0 + R_0)**k)/(4.0*(1.0 + R)**k))

def plot_profile(z, values, x_label, y_label, filename, xticks=[], yticks=[], y_inverted = False, y_log = False, y_lim = []):
    #xticks = [x_start, x_end, x_num]
    plt.rc('text', usetex=latex_labels)

    fig = plt.figure()

    plt.plot(values, z)
    plt.xlabel(x_label)
    plt.ylabel(y_label)

    if y_log:
        plt.yscale('log')
    if y_lim:
        plt.gca().set_ylim([y_lim[0],y_lim[1]])
    if y_inverted:
        plt.gca().invert_yaxis()
    if xticks:
        plt.xticks(np.linspace(xticks[0],xticks[1],xticks[2],endpoint=True))
        plt.gca().get_xaxis().get_major_formatter().labelOnlyBase = False
        plt.gca().get_xaxis().set_major_formatter(mpl.ticker.FormatStrFormatter('%d'))
    if yticks:
        plt.yticks(np.linspace(yticks[0],yticks[1],yticks[2],endpoint=True))
        plt.gca().get_yaxis().get_major_formatter().labelOnlyBase = False
        plt.gca().get_yaxis().set_major_formatter(mpl.ticker.FormatStrFormatter('%d'))
    plt.savefig(filename, bbox_inches = 'tight', pad_inches = 0.1)
    plt.close()

def plot_profile_multi(z, values, labels, x_label, y_label, filename, obs_z=None, obs_values=None, xticks=[], yticks=[], y_inverted = False, y_log = False, y_lim = [], conversion_factor=1.0, freeze_axis=None, line_type=None, color_index=None, skill_scores=False, zero_line=False):
    #xticks = [x_start, x_end, x_num]
    plt.rc('text', usetex=latex_labels)
    
    if np.count_nonzero(values) == 0:
        print('The plot named {} will not be created due to all zero values'.format(filename))
        return
    
    fig = plt.figure()
    colors = ['#e41a1c','#4daf4a','#377eb8','#984ea3','#ff7f00','#a65628','#f781bf','#ffff33']
    linestyles = ['-','--','-.',':']
    markers = ['+','o','*','s','d','^','v','p','h']
    linewidth_val = 2.0
    multi_legend = True
    if multi_legend and freeze_axis is None:
        legend1_labels = []
        legend2_labels = []
        lines1 = []
        lines2 = []
    else:
        legend_labels = []
        lines = []

    line_types = ['color', 'style']
    if line_type is None:
        line_type = 'color'
    else:
        if line_type not in line_types:
            raise ValueError("Invalid line type. Expected one of: %s" % line_types)
    if color_index is None:
        color_index = 0
    else:
        if color_index > len(colors):
            raise ValueError("Invalid color index. Expected color index < " + string(len(colors)))
        else:
            #move color_index color to the front of the colors list (leave rest of list in same order)
            colors.insert(0, colors.pop(color_index))


    if any(isinstance(i, list) for i in values): #check for a list of lists
        if freeze_axis is None:
            if multi_legend:
                # for i in range(len(values)):
                #     lines1.append(plt.Line2D((0,1),(0,0), linestyle=linestyles[0], color=colors[i], linewidth=linewidth_val))
                #     legend1_labels.append(labels[0][i])
                # for j in range(len(values[0])):
                #     if j < len(linestyles):
                #         lines2.append(plt.Line2D((0,1),(0,0), linestyle=linestyles[j], color='black', linewidth=linewidth_val))
                #     else:
                #         lines2.append(plt.Line2D((0,1),(0,0), linestyle='', marker=markers[j - len(linestyles)], color='black', linewidth=linewidth_val))
                #     legend2_labels.append(labels[1][j])
                for i in range(len(values)):
                    if i < len(linestyles):
                        lines1.append(plt.Line2D((0,1),(0,0), linestyle=linestyles[i], color='black', linewidth=linewidth_val))
                    else:
                        lines1.append(plt.Line2D((0,1),(0,0), linestyle='', marker=markers[i - len(linestyles)], color='black', linewidth=linewidth_val))
                    legend1_labels.append(labels[0][i])
                for j in range(len(values[0])):
                    lines2.append(plt.Line2D((0,1),(0,0), color=colors[j], linewidth=linewidth_val))
                    legend2_labels.append(labels[1][j])

            # for i in range(len(values)): #number of lists in the list (linestyles)
            #     for j in range(len(values[i])): #number of items in the nested lists (colors)
            #         if j < len(linestyles):
            #             plt.plot(values[i][j], z, linestyles[j], color=colors[i], linewidth=linewidth_val)
            #         else:
            #             plt.plot(values[i][j], z, markers[j - len(linestyles)], color=colors[i], linewidth=linewidth_val)
            #         if not multi_legend:
            #             if j < len(linestyles):
            #                 lines.append(plt.Line2D((0,1),(0,0),color=colors[i],linestyle=linestyles[j], linewidth=linewidth_val))
            #             else:
            #                 lines.append(plt.Line2D((0,1),(0,0),linestyle='', color=colors[i], marker=markers[j - len(linestyles)], linewidth=linewidth_val))
            #             legend_labels.append(labels[0][i] + '-' + labels[1][j])
            for i in range(len(values)): #number of lists in the list (linestyles)
                for j in range(len(values[i])): #number of items in the nested lists (colors)
                    if i < len(linestyles):
                        plt.plot(conversion_factor*values[i][j], z, linestyles[i], color=colors[j], linewidth=linewidth_val)
                    else:
                        plt.plot(conversion_factor*values[i][j], z, markers[i - len(linestyles)], color=colors[j], linewidth=linewidth_val)
                    if not multi_legend:
                        if i < len(linestyles):
                            lines.append(plt.Line2D((0,1),(0,0),color=colors[j],linestyle=linestyles[i], linewidth=linewidth_val))
                        else:
                            lines.append(plt.Line2D((0,1),(0,0),linestyle='', color=colors[j],marker=markers[i - len(linestyles)], linewidth=linewidth_val))
                        legend_labels.append(labels[0][i] + '-' + labels[1][j])

        else:
            legend_labels = []
            for i in range(len(values)): #number of lists in the list (colors)
                plt.plot(conversion_factor*values[i][freeze_axis], z, color=colors[i], linestyle=linestyles[0], linewidth=linewidth_val)
                lines.append(plt.Line2D((0,1),(0,0),color=colors[i],linestyle=linestyles[0], linewidth=linewidth_val))
                legend_labels.append(labels[0][i])
    else:
        multi_legend = False
        legend_labels = labels
        lines = []
        ss = []
        for i in range(len(values)):
            if line_type == 'style':
                if i < len(linestyles):
                    plt.plot(conversion_factor*values[i], z, color=colors[0], linestyle=linestyles[i], linewidth=linewidth_val)
                    lines.append(plt.Line2D((0,1),(0,0), color=colors[0], linestyle=linestyles[i], linewidth=linewidth_val))
                else:
                    plt.plot(conversion_factor*values[i], z, color=colors[0], linestyle='', marker=markers[i - len(linestyles)], linewidth=linewidth_val)
                    lines.append(plt.Line2D((0,1),(0,0),color=colors[0], linestyle='', marker=markers[i - len(linestyles)], linewidth=linewidth_val))
            else:
                plt.plot(conversion_factor*values[i], z, color=colors[i], linestyle=linestyles[0], linewidth=linewidth_val)
                lines.append(plt.Line2D((0,1),(0,0), color=colors[i], linestyle=linestyles[0], linewidth=linewidth_val))
            if skill_scores and obs_values is not None and obs_z is not None:
                interp_values = np.flipud(np.interp(np.flipud(obs_z), np.flipud(z), np.flipud(values[i])))
                #print interp_values, obs_z, z
                ss.append(calc_skill_score(interp_values, obs_values, 1.0, 4.0))
                #print ss
                #legend_labels[i] = legend_labels[i] + '(' + str(ss) + ')'

    if obs_values is not None and obs_z is not None:
        plt.plot(conversion_factor*obs_values, obs_z, color='black', linewidth=linewidth_val)
        lines.append(plt.Line2D((0,1),(0,0),color='black', linewidth=linewidth_val))
        legend_labels.append('obs')

    if zero_line:
        plt.axvline(x=0.0, zorder=1, color='black')

    plt.xlabel(x_label)
    plt.ylabel(y_label)

    if multi_legend and freeze_axis is None:
        first_legend = plt.legend(lines2, legend2_labels, loc='upper left', fontsize='small', bbox_to_anchor=(1.01, 1), borderaxespad=0.)
        ax = plt.gca().add_artist(first_legend)
        plt.legend(lines1, legend1_labels, loc='lower left', fontsize='small', bbox_to_anchor=(1.01, 0), borderaxespad=0.)
    else:
        if skill_scores and obs_values is not None and obs_z is not None:
            updated_legend_labels = []
            for i in range(len(values)):
                updated_legend_labels.append(legend_labels[i] + ' (' + '{0:.3f}'.format(ss[i]) + ')')
            first_legend = plt.legend(lines, updated_legend_labels, loc='best', fontsize='small')
        else:
            first_legend = plt.legend(lines, legend_labels, loc='best', fontsize='small')

    if y_log:
        plt.yscale('log')
    if y_lim:
        plt.gca().set_ylim([y_lim[0],y_lim[1]])
        #change x axis to suit reduced y axis
        min_x = []
        max_x = []
        for i in range(len(values)):
            include_y = np.intersect1d(np.where(z >= y_lim[0])[0], np.where(z <= y_lim[1])[0])
            start_y = include_y[0]
            end_y = include_y[-1]
            min_x.append(conversion_factor*np.min(values[i][start_y:end_y]))
            max_x.append(conversion_factor*np.max(values[i][start_y:end_y]))
        if obs_values is not None and obs_z is not None:
            include_y = np.intersect1d(np.where(obs_z >= y_lim[0])[0], np.where(obs_z <= y_lim[1])[0])
            start_y = include_y[0]
            end_y = include_y[-1]
            min_x.append(conversion_factor*np.min(obs_values[start_y:end_y]))
            max_x.append(conversion_factor*np.max(obs_values[start_y:end_y]))
        min_x_all = min(min_x)
        max_x_all = max(max_x)
        plt.gca().set_xlim([min_x_all, max_x_all])
    if y_inverted:
        plt.gca().invert_yaxis()
    if xticks:
        plt.xticks(np.linspace(xticks[0],xticks[1],xticks[2],endpoint=True))
        plt.gca().get_xaxis().get_major_formatter().labelOnlyBase = False
        plt.gca().get_xaxis().set_major_formatter(mpl.ticker.FormatStrFormatter('%d'))
    if yticks:
        plt.yticks(np.linspace(yticks[0],yticks[1],yticks[2],endpoint=True))
        plt.gca().get_yaxis().get_major_formatter().labelOnlyBase = False
        plt.gca().get_yaxis().set_major_formatter(mpl.ticker.FormatStrFormatter('%d'))
    plt.savefig(filename, bbox_inches = 'tight')
    plt.close()

def plot_profile_compare(z, values, LES_values, LES_z, x_label, y_label, filename, xticks=[], yticks=[], y_inverted = False, y_log = False, x_lim = [], y_lim = []):
    #xticks = [x_start, x_end, x_num]
    plt.rc('text', usetex=latex_labels)

    fig = plt.figure()

    SCM_color = 'red'
    LES_color = 'black'

    plt.plot(values, z, color=SCM_color, linewidth=2)
    plt.plot(LES_values[0,:], LES_z, color=LES_color, linestyle=':')
    plt.plot(LES_values[1,:], LES_z, color=LES_color, linestyle='--')
    plt.plot(LES_values[2,:], LES_z, color=LES_color, linestyle=':')
    plt.xlabel(x_label)
    plt.ylabel(y_label)

    if y_log:
        plt.yscale('log')
    if x_lim:
        plt.gca().set_xlim([x_lim[0],x_lim[1]])
    if y_lim:
        plt.gca().set_ylim([y_lim[0],y_lim[1]])
    if y_inverted:
        plt.gca().invert_yaxis()
    if xticks:
        plt.xticks(np.linspace(xticks[0],xticks[1],xticks[2],endpoint=True))
        plt.gca().get_xaxis().get_major_formatter().labelOnlyBase = False
        plt.gca().get_xaxis().set_major_formatter(mpl.ticker.FormatStrFormatter('%d'))
    if yticks:
        plt.yticks(np.linspace(yticks[0],yticks[1],yticks[2],endpoint=True))
        plt.gca().get_yaxis().get_major_formatter().labelOnlyBase = False
        plt.gca().get_yaxis().set_major_formatter(mpl.ticker.FormatStrFormatter('%d'))
    plt.savefig(filename, bbox_inches = 'tight', pad_inches = 0.1)
    plt.close()

def plot_time_series(time, values, x_label, y_label, filename):
    plt.rc('text', usetex=latex_labels)

    fig = plt.figure()

    plt.plot(time, values)
    plt.xlabel(x_label)
    plt.ylabel(y_label)

    plt.savefig(filename, bbox_inches = 'tight', pad_inches = 0.1)
    plt.close()

def plot_time_series_compare(time, values, LES_time, LES_values, x_label, y_label, filename):
    plt.rc('text', usetex=latex_labels)

    fig = plt.figure()

    SCM_color = 'red'
    LES_color = 'black'

    plt.plot(time, values)
    plt.plot(LES_time, LES_values[0,:], color=LES_color, linestyle=':')
    plt.plot(LES_time, LES_values[1,:], color=LES_color, linestyle='--')
    plt.plot(LES_time, LES_values[2,:], color=LES_color, linestyle=':')

    plt.xlabel(x_label)
    plt.ylabel(y_label)

    plt.savefig(filename, bbox_inches = 'tight', pad_inches = 0.1)
    plt.close()

def plot_time_series_multi(time, values, labels, x_label, y_label, filename, obs_time=None, obs_values=None, obs_label=None, line_type=None, color_index=None, skill_scores=False, conversion_factor=1.0):
    plt.rc('text', usetex=latex_labels)
    
    if np.count_nonzero(values) == 0:
        print('The plot named {} will not be created due to all zero values'.format(filename))
        return
    
    fig = plt.figure()
    
    lines = []
    colors = ['#e41a1c','#4daf4a','#377eb8','#984ea3','#ff7f00','#a65628','#f781bf','#ffff33']
    linestyles = ['-','--','-.',':']
    markers = ['+','o','*','s','d','^','v','p','h']
    linewidth_val = 2.0
    multi_legend = True
    if multi_legend:
        legend1_labels = []
        legend2_labels = []
        lines1 = []
        lines2 = []
    else:
        legend_labels = []
        lines = []

    line_types = ['color', 'style']
    if line_type is None:
        line_type = 'color'
    else:
        if line_type not in line_types:
            raise ValueError("Invalid line type. Expected one of: %s" % line_types)

    if color_index is None:
        color_index = 0
    else:
        if color_index > len(colors):
            raise ValueError("Invalid color index. Expected color index < " + string(len(colors)))
        else:
            #move color_index color to the front of the colors list (leave rest of list in same order)
            colors.insert(0, colors.pop(color_index))

    if any(isinstance(i, list) for i in values): #check for a list of lists
        if multi_legend:
            for i in range(len(values)):
                if i < len(linestyles):
                    lines1.append(plt.Line2D((0,1),(0,0),color='black',linestyle=linestyles[i], linewidth=linewidth_val))
                else:
                    lines1.append(plt.Line2D((0,1),(0,0),color='black',linestyle='', marker=markers[i - len(linestyles)], linewidth=linewidth_val))
                legend1_labels.append(labels[0][i])
            for j in range(len(values[0])):
                lines2.append(plt.Line2D((0,1),(0,0),color=colors[j],linestyle='-', linewidth=linewidth_val))
                legend2_labels.append(labels[1][j])

        for i in range(len(values)): #number of lists in the list (linestyles/markers)
            for j in range(len(values[i])): #number of items in the nested lists (colors)
                if i < len(linestyles):
                    plt.plot(time, conversion_factor*values[i][j], color=colors[j], linestyle=linestyles[i], linewidth=linewidth_val)
                    lines.append(plt.Line2D((0,1),(0,0),color=colors[j],linestyle=linestyles[i], linewidth=linewidth_val))
                else:
                    plt.plot(time, conversion_factor*values[i][j], color=colors[j], linestyle='', marker=markers[i - len(linestyles)], linewidth=linewidth_val)
                    lines.append(plt.Line2D((0,1),(0,0),color=colors[j],linestyle='', marker=markers[i - len(linestyles)], linewidth=linewidth_val))
    else:
        multi_legend = False
        legend_labels = labels
        lines = []
        ss = []
        for i in range(len(values)):
            if line_type == 'style':
                if i < len(linestyles):
                    plt.plot(time, conversion_factor*values[i], color=colors[0], linestyle=linestyles[i], linewidth=linewidth_val)
                    lines.append(plt.Line2D((0,1),(0,0),color=colors[0], linestyle=linestyles[i], linewidth=linewidth_val))
                else:
                    plt.plot(time, conversion_factor*values[i], color=colors[0], linestyle='', marker = markers[i - len(linestyles)], linewidth=linewidth_val)
                    lines.append(plt.Line2D((0,1),(0,0),color=colors[0], linestyle='', marker = markers[i - len(linestyles)], linewidth=linewidth_val))
            else:
                plt.plot(time, conversion_factor*values[i], color=colors[i], linestyle=linestyles[0], linewidth=linewidth_val)
                lines.append(plt.Line2D((0,1),(0,0),color=colors[i], linestyle=linestyles[0], linewidth=linewidth_val))
            if skill_scores and obs_values is not None and obs_time is not None:
                #print obs_time, time
                #interp_values = np.interp(obs_time, time, values[i])
                #print interp_values, obs_z, z
                ss.append(calc_skill_score(values[i], obs_values, 1.0, 4.0))
                #print ss
                #legend_labels[i] = legend_labels[i] + '(' + str(ss) + ')'
    if obs_values is not None and obs_time is not None:
        plt.plot(obs_time, conversion_factor*obs_values, color='black', linewidth=linewidth_val)
        lines.append(plt.Line2D((0,1),(0,0),color='black', linewidth=linewidth_val))
        if obs_label is not None:
            legend_labels.append(obs_label)
        else:
            legend_labels.append('obs')

    #plt.legend(lines, labels, loc='best', fontsize='small')

    plt.xlabel(x_label)
    plt.ylabel(y_label)

    if multi_legend:
        first_legend = plt.legend(lines2, legend2_labels, loc='upper left', fontsize='small', bbox_to_anchor=(1.01, 1), borderaxespad=0.)
        ax = plt.gca().add_artist(first_legend)
        plt.legend(lines1, legend1_labels, loc='lower left', fontsize='small', bbox_to_anchor=(1.01, 0), borderaxespad=0.)
    else:
        if skill_scores and obs_values is not None and obs_time is not None:
            updated_legend_labels = []
            for i in range(len(values)):
                updated_legend_labels.append(legend_labels[i] + ' (' + '{0:.3f}'.format(ss[i]) + ')')
            first_legend = plt.legend(lines, updated_legend_labels, loc='best', fontsize='small')
        else:
            first_legend = plt.legend(lines, legend_labels, loc='best', fontsize='small')

    plt.savefig(filename, bbox_inches = 'tight', pad_inches = 0.1)
    plt.close()

def contour_plot_firl(x_dim, y_dim, values, min_val, max_val, title, x_label, y_label, filename, xticks=[], yticks=[], plot_mean = 0, annotation = 0, y_inverted = 0, y_log = False, y_lim = [], conversion_factor=1.0):
    #os.environ['PATH'] = os.environ['PATH'] + ':/usr/texbin'
    #plt.rc('ps', usedistiller='xpdf')
    plt.rc('text', usetex=latex_labels)
    
    if np.count_nonzero(values) == 0:
        print('The plot named {} will not be created due to all zero values'.format(filename))
        return
        
    if np.amax(values) == np.amin(values):
        print('The plot named {} will not be created due to all values being equal'.format(filename))
        return
    
    if(min_val != -999 and max_val != -999):
        min_val = conversion_factor*min_val
        max_val = conversion_factor*max_val
        if min_val < 0.0 and max_val > 0.0:
            colormap = 'bwr'
            min_contour_value = min_val
            max_contour_value = max_val
        else:
            colormap = 'gist_yarg'
            min_contour_value = min_val
            max_contour_value = max_val
    else:
        if np.min(values) < 0.0 and np.max(values) > 0.0:
            colormap = 'bwr'
            max_abs = np.max(np.abs(conversion_factor*values))
            min_contour_value = -1.0*max_abs
            max_contour_value = max_abs
        else:
            colormap = 'gist_yarg'
            min_contour_value = np.min(conversion_factor*values)
            max_contour_value = np.max(conversion_factor*values)

    fig = plt.figure()

    #plt.contourf(x_dim, y_dim, values, 8, alpha=.75, cmap='jet')
    v = np.linspace(min_contour_value, max_contour_value, 8, endpoint=True)
    plt.contourf(x_dim, y_dim, conversion_factor*values, v, vmin=min_contour_value, vmax = max_contour_value, alpha=1.0, cmap=colormap)
    cb = plt.colorbar(ticks=v)
    C = plt.contour(x_dim, y_dim, conversion_factor*values, v, vmin=min_contour_value, vmax = max_contour_value, colors='black')
    #contour_labels = plt.clabel(C, inline=1, fontsize=10) #when saving as eps, contour labels are outside of the axes bounding box

    if annotation:
        ax = plt.gca()
        ax.text(0.6*x_end, 0.9*y_end, annotation)

    #plt.xticks(np.linspace(x_start,x_end,x_num,endpoint=True)), plt.yticks(np.linspace(y_start,y_end,y_num,endpoint=True))

    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)

    #plt.tight_layout()

    if plot_mean:
        plt.subplots_adjust(left = 0.1, bottom = 0.1, right=0.85, top = 0.9, wspace=0.0)
        mean_ax = fig.add_axes([0.85, 0.10, 0.13, 0.8])
        mean_profile = np.mean(conversion_factor*values, axis=1, dtype=np.float64)
        plt.plot(mean_profile, y_dim)
        plt.xticks(np.linspace(np.min(mean_profile),np.max(mean_profile),2,endpoint=True)), plt.yticks(np.linspace(y_start,y_end,y_num,endpoint=True))
        plt.setp( mean_ax.get_yticklabels(), visible=False)
        if(min_val != -999 and max_val != -999):
            ax = plt.gca()
            ax.set_xlim((min_val,max_val))

    if y_log:
        plt.yscale('log')
    if y_lim:
        plt.gca().set_ylim([y_lim[0],y_lim[1]])
    if y_inverted:
        plt.gca().invert_yaxis()
    if xticks:
        plt.xticks(np.linspace(xticks[0],xticks[1],xticks[2],endpoint=True))
        plt.gca().get_xaxis().get_major_formatter().labelOnlyBase = False
        plt.gca().get_xaxis().set_major_formatter(mpl.ticker.FormatStrFormatter('%d'))
    if yticks:
        plt.yticks(np.linspace(yticks[0],yticks[1],yticks[2],endpoint=True))
        plt.gca().get_yaxis().get_major_formatter().labelOnlyBase = False
        plt.gca().get_yaxis().set_major_formatter(mpl.ticker.FormatStrFormatter('%d'))

    plt.savefig(filename, bbox_inches = 'tight', pad_inches = 0.1)
    plt.close()

def plot_profile_multi_ens(z, values, labels, x_label, y_label, filename, obs_z=None, obs_values=None, xticks=[], yticks=[], y_inverted = False, y_log = False, y_lim = [], freeze_axis=None, line_type=None, color_index=None, skill_scores=False):
    #xticks = [x_start, x_end, x_num]
    plt.rc('text', usetex=latex_labels)

    #ensemble data processing
    values_ens_processed = []
    if any(isinstance(i, list) for i in values): #check for a list of lists
        for i in range(len(values)):
            dummy_list = []
            for j in range(len(values[i])):
                dummy_list.append(np.percentile(values[i][j], [0, 10, 25, 50, 75, 90, 100], axis=0, interpolation='nearest'))
            values_ens_processed.append(dummy_list)
    else:
        for i in range(len(values)):
            values_ens_processed.append(np.percentile(values[i], [0, 10, 25, 50, 75, 90, 100], axis=0, interpolation='nearest'))


    fig = plt.figure()
    colors = ['#e41a1c','#4daf4a','#377eb8','#984ea3','#ff7f00','#a65628','#f781bf','#ffff33']
    linestyles = ['-','--','-.',':']
    markers = ['+','o','*','s','d','^','v','p','h']
    hatches = ['/','\\','|','-']
    linewidth_val = 3.0

    multi_legend = True
    if multi_legend and freeze_axis is None:
        legend1_labels = []
        legend2_labels = []
        lines1 = []
        lines2 = []
    else:
        legend_labels = []
        lines = []

    line_types = ['color', 'style']
    if line_type is None:
        line_type = 'color'
    else:
        if line_type not in line_types:
            raise ValueError("Invalid line type. Expected one of: %s" % line_types)

    if color_index is None:
        color_index = 0
    else:
        if color_index > len(colors):
            raise ValueError("Invalid color index. Expected color index < " + string(len(colors)))
        else:
            #move color_index color to the front of the colors list (leave rest of list in same order)
            colors.insert(0, colors.pop(color_index))


    if any(isinstance(i, list) for i in values): #check for a list of lists
        if freeze_axis is None:
            if multi_legend:
                for i in range(len(values)):
                    if i < len(linestyles):
                        lines1.append(plt.Line2D((0,1),(0,0), linestyle=linestyles[i], color='black', linewidth=linewidth_val))
                    else:
                        lines1.append(plt.Line2D((0,1),(0,0), linestyle='', marker=markers[i - len(linestyles)], color='black', linewidth=linewidth_val))
                    legend1_labels.append(labels[0][i])
                for j in range(len(values[0])):
                    lines2.append(plt.Line2D((0,1),(0,0), color=colors[j], linewidth=linewidth_val))
                    legend2_labels.append(labels[1][j])

            for i in range(len(values)): #number of lists in the list (linestyles)
                for j in range(len(values[i])): #number of items in the nested lists (colors)
                    if i < len(linestyles):
                        #first plot the range (10-90 percentiles)
                        plt.fill_betweenx(z, values_ens_processed[i][j][1,:], values_ens_processed[i][j][-2,:], facecolor=colors[j], alpha=0.25, hatch=hatches[i])
                        #plot the median
                        plt.plot(values_ens_processed[i][j][3,:], z, color=colors[j], linestyle=linestyles[i], linewidth=linewidth_val)
                    else:
                        #first plot the range (10-90 percentiles)
                        plt.fill_betweenx(z, values_ens_processed[i][j][1,:], values_ens_processed[i][j][-2,:], facecolor=colors[j], alpha=0.25, hatch=markers[i - len(linestyles)])
                        #plot the median
                        plt.plot(values_ens_processed[i][j][3,:], z, color=colors[j], linestyle='', marker=markers[i - len(linestyles)], linewidth=linewidth_val)

                    if not multi_legend:
                        if i < len(linestyles):
                            lines.append(plt.Line2D((0,1),(0,0),color=colors[j],linestyle=linestyles[i], linewidth=linewidth_val))
                        else:
                            lines.append(plt.Line2D((0,1),(0,0),linestyle='', color=colors[j],marker=markers[i - len(linestyles)], linewidth=linewidth_val))
                        legend_labels.append(labels[0][i] + '-' + labels[1][j])

        else:
            legend_labels = []
            for i in range(len(values)): #number of lists in the list (colors)
                plt.plot(values[i][freeze_axis], z, color=colors[i], linestyle=linestyles[0], linewidth=linewidth_val)
                lines.append(plt.Line2D((0,1),(0,0),color=colors[i],linestyle=linestyles[0], linewidth=linewidth_val))
                legend_labels.append(labels[0][i])
    else:
        multi_legend = False
        legend_labels = labels
        lines = []
        ss = []
        for i in range(len(values)):
            if line_type == 'style':
                if i < len(linestyles):
                    #first plot the range (10-90 percentiles)
                    plt.fill_betweenx(z, values_ens_processed[i][1,:], values_ens_processed[i][-2,:], facecolor=colors[0], alpha=0.25, hatch=hatches[i])
                    #plot the 5th and 95th percentiles
                    # plt.plot(values_ens_processed[i][1,:], z, color=colors[0], linestyle=linestyles[i], linewidth=0.05*linewidth_val)
                    # plt.plot(values_ens_processed[i][-2,:], z, color=colors[0], linestyle=linestyles[i], linewidth=0.05*linewidth_val)
                    #plot the 25th and 75th percentiles
                    # plt.plot(values_ens_processed[i][2,:], z, color=colors[0], linestyle=linestyles[i], linewidth=0.25*linewidth_val)
                    # plt.plot(values_ens_processed[i][-3,:], z, color=colors[0], linestyle=linestyles[i], linewidth=0.25*linewidth_val)
                    #plot the median
                    plt.plot(values_ens_processed[i][3,:], z, color=colors[0], linestyle=linestyles[i], linewidth=linewidth_val)
                    lines.append(plt.Line2D((0,1),(0,0), color=colors[0], linestyle=linestyles[i], linewidth=linewidth_val))
                else:
                    #first plot the range (10-90 percentiles)
                    plt.fill_betweenx(z, values_ens_processed[i][1,:], values_ens_processed[i][-2,:], facecolor=colors[0], alpha=0.25, hatch=markers[i - len(linestyles)])
                    #plot the 5th and 95th percentiles
                    # plt.plot(values_ens_processed[i][1,:], z, color=colors[0], linestyle='', marker=markers[i - len(linestyles)], linewidth=0.05*linewidth_val)
                    # plt.plot(values_ens_processed[i][-2,:], z, color=colors[0], linestyle='', marker=markers[i - len(linestyles)], linewidth=0.05*linewidth_val)
                    #plot the 25th and 75th percentiles
                    # plt.plot(values_ens_processed[i][2,:], z, color=colors[0], linestyle='', marker=markers[i - len(linestyles)], linewidth=0.25*linewidth_val)
                    # plt.plot(values_ens_processed[i][-3,:], z, color=colors[0], linestyle='', marker=markers[i - len(linestyles)], linewidth=0.25*linewidth_val)
                    #plot the median
                    plt.plot(values_ens_processed[i][3,:], z, color=colors[0], linestyle='', marker=markers[i - len(linestyles)], linewidth=linewidth_val)
                    lines.append(plt.Line2D((0,1),(0,0),color=colors[0], linestyle='', marker=markers[i - len(linestyles)], linewidth=linewidth_val))
            else:
                #first plot the range (10-90 percentiles)
                plt.fill_betweenx(z, values_ens_processed[i][1,:], values_ens_processed[i][-2,:], facecolor=colors[i], alpha=0.25)
                #plot the 5th and 95th percentiles
                # plt.plot(values_ens_processed[i][1,:], z, color=colors[i], linestyle=linestyles[0], linewidth=0.05*linewidth_val)
                # plt.plot(values_ens_processed[i][-2,:], z, color=colors[i], linestyle=linestyles[0], linewidth=0.05*linewidth_val)
                #plot the 25th and 75th percentiles
                plt.plot(values_ens_processed[i][2,:], z, color=colors[i], linestyle=linestyles[0], linewidth=0.25*linewidth_val)
                plt.plot(values_ens_processed[i][-3,:], z, color=colors[i], linestyle=linestyles[0], linewidth=0.25*linewidth_val)
                #plot the median
                plt.plot(values_ens_processed[i][3,:], z, color=colors[i], linestyle=linestyles[0], linewidth=linewidth_val)
                lines.append(plt.Line2D((0,1),(0,0), color=colors[i], linestyle=linestyles[0], linewidth=linewidth_val))
            if skill_scores and obs_values is not None and obs_z is not None:
                interp_values = np.flipud(np.interp(np.flipud(obs_z), np.flipud(z), np.flipud(values_ens_processed[i][3,:])))
                #print interp_values, obs_z, z
                ss.append(calc_skill_score(interp_values, obs_values, 1.0, 4.0))
                #print ss
                #legend_labels[i] = legend_labels[i] + '(' + str(ss) + ')'
    if obs_values is not None and obs_z is not None:
        plt.plot(obs_values, obs_z, color='black', linewidth=linewidth_val)
        lines.append(plt.Line2D((0,1),(0,0),color='black', linewidth=linewidth_val))
        legend_labels.append('obs')

    plt.xlabel(x_label)
    plt.ylabel(y_label)

    if multi_legend and freeze_axis is None:
        first_legend = plt.legend(lines2, legend2_labels, loc='upper left', fontsize='small', bbox_to_anchor=(1.01, 1), borderaxespad=0.)
        ax = plt.gca().add_artist(first_legend)
        plt.legend(lines1, legend1_labels, loc='lower left', fontsize='small', bbox_to_anchor=(1.01, 0), borderaxespad=0.)
    else:
        if skill_scores and obs_values is not None and obs_z is not None:
            updated_legend_labels = []
            for i in range(len(values)):
                updated_legend_labels.append(legend_labels[i] + ' (' + '{0:.3f}'.format(ss[i]) + ')')
            first_legend = plt.legend(lines, updated_legend_labels, loc='best', fontsize='small')
        else:
            first_legend = plt.legend(lines, legend_labels, loc='best', fontsize='small')

    if y_log:
        plt.yscale('log')
    if y_lim:
        plt.gca().set_ylim([y_lim[0],y_lim[1]])
    if y_inverted:
        plt.gca().invert_yaxis()
    if xticks:
        plt.xticks(np.linspace(xticks[0],xticks[1],xticks[2],endpoint=True))
        plt.gca().get_xaxis().get_major_formatter().labelOnlyBase = False
        plt.gca().get_xaxis().set_major_formatter(mpl.ticker.FormatStrFormatter('%d'))
    if yticks:
        plt.yticks(np.linspace(yticks[0],yticks[1],yticks[2],endpoint=True))
        plt.gca().get_yaxis().get_major_formatter().labelOnlyBase = False
        plt.gca().get_yaxis().set_major_formatter(mpl.ticker.FormatStrFormatter('%d'))
    plt.savefig(filename, bbox_inches = 'tight')
    plt.close()

def plot_time_series_multi_ens(time, values, labels, x_label, y_label, filename, obs_time=None, obs_values=None, obs_label=None, line_type=None, color_index=None, skill_scores=False):
    plt.rc('text', usetex=latex_labels)

    fig = plt.figure()

    #ensemble data processing
    values_ens_processed = []
    if any(isinstance(i, list) for i in values): #check for a list of lists
        for i in range(len(values)):
            dummy_list = []
            for j in range(len(values[i])):
                dummy_list.append(np.percentile(values[i][j], [0, 10, 25, 50, 75, 90, 100], axis=0, interpolation='nearest'))
            values_ens_processed.append(dummy_list)
    else:
        for i in range(len(values)):
            values_ens_processed.append(np.percentile(values[i], [0, 10, 25, 50, 75, 90, 100], axis=0, interpolation='nearest'))

    lines = []
    colors = ['#e41a1c','#4daf4a','#377eb8','#984ea3','#ff7f00','#a65628','#f781bf','#ffff33']
    linestyles = ['-','--','-.',':']
    markers = ['+','o','*','s','d','^','v','p','h']
    hatches = ['/','\\','|','-']
    linewidth_val = 2.0
    multi_legend = True
    if multi_legend:
        legend1_labels = []
        legend2_labels = []
        lines1 = []
        lines2 = []
    else:
        legend_labels = []
        lines = []

    line_types = ['color', 'style']
    if line_type is None:
        line_type = 'color'
    else:
        if line_type not in line_types:
            raise ValueError("Invalid line type. Expected one of: %s" % line_types)

    if color_index is None:
        color_index = 0
    else:
        if color_index > len(colors):
            raise ValueError("Invalid color index. Expected color index < " + string(len(colors)))
        else:
            #move color_index color to the front of the colors list (leave rest of list in same order)
            colors.insert(0, colors.pop(color_index))

    if any(isinstance(i, list) for i in values): #check for a list of lists
        if multi_legend:
            for i in range(len(values)):
                if i < len(linestyles):
                    lines1.append(plt.Line2D((0,1),(0,0),color='black',linestyle=linestyles[i], linewidth=linewidth_val))
                else:
                    lines1.append(plt.Line2D((0,1),(0,0),color='black',linestyle='', marker=markers[i - len(linestyles)], linewidth=linewidth_val))
                legend1_labels.append(labels[0][i])
            for j in range(len(values[0])):
                lines2.append(plt.Line2D((0,1),(0,0),color=colors[j],linestyle='-', linewidth=linewidth_val))
                legend2_labels.append(labels[1][j])

        for i in range(len(values)): #number of lists in the list (linestyles/markers)
            for j in range(len(values[i])): #number of items in the nested lists (colors)
                if i < len(linestyles):
                    #first plot the range (10-90 percentiles)
                    plt.fill_between(time, values_ens_processed[i][j][1,:], values_ens_processed[i][j][-2,:], facecolor=colors[j], alpha=0.25, hatch=hatches[i])
                    #plot the median
                    plt.plot(time, values_ens_processed[i][j][3,:], color=colors[j], linestyle=linestyles[i], linewidth=linewidth_val)

                    #plt.plot(time, values[i][j], color=colors[j], linestyle=linestyles[i], linewidth=linewidth_val)
                    lines.append(plt.Line2D((0,1),(0,0),color=colors[j],linestyle=linestyles[i], linewidth=linewidth_val))
                else:
                    #first plot the range (10-90 percentiles)
                    plt.fill_between(time, values_ens_processed[i][j][1,:], values_ens_processed[i][j][-2,:], facecolor=colors[j], alpha=0.25, hatch=markers[i - len(linestyles)])
                    #plot the median
                    plt.plot(time, values_ens_processed[i][j][3,:], color=colors[j], linestyle='', marker=markers[i - len(linestyles)], linewidth=linewidth_val)

                    #plt.plot(time, values[i][j], color=colors[j], linestyle='', marker=markers[i - len(linestyles)], linewidth=linewidth_val)
                    lines.append(plt.Line2D((0,1),(0,0),color=colors[j],linestyle='', marker=markers[i - len(linestyles)], linewidth=linewidth_val))
    else:
        multi_legend = False
        legend_labels = labels
        lines = []
        ss = []
        for i in range(len(values)):
            if line_type == 'style':
                if i < len(linestyles):
                    #first plot the range (10-90 percentiles)
                    plt.fill_between(time, values_ens_processed[i][1,:], values_ens_processed[i][-2,:], facecolor=colors[0], alpha=0.25, hatch=hatches[i])
                    #plot the median
                    plt.plot(time, values_ens_processed[i][3,:], color=colors[0], linestyle=linestyles[i], linewidth=linewidth_val)

                    lines.append(plt.Line2D((0,1),(0,0),color=colors[0], linestyle=linestyles[i], linewidth=linewidth_val))
                else:
                    #first plot the range (10-90 percentiles)
                    plt.fill_between(time, values_ens_processed[i][1,:], values_ens_processed[i][-2,:], facecolor=colors[0], alpha=0.25, hatch=markers[i - len(linestyles)])
                    #plot the median
                    plt.plot(time, values_ens_processed[i][3,:], color=colors[0], linestyle='', marker=markers[i - len(linestyles)], linewidth=linewidth_val)

                    lines.append(plt.Line2D((0,1),(0,0),color=colors[0], linestyle='', marker = markers[i - len(linestyles)], linewidth=linewidth_val))
            else:
                #first plot the range (10-90 percentiles)
                plt.fill_between(time, values_ens_processed[i][1,:], values_ens_processed[i][-2,:], facecolor=colors[i], alpha=0.25)
                #plot the 25th and 75th percentiles
                plt.plot(time, values_ens_processed[i][2,:], color=colors[i], linestyle=linestyles[0], linewidth=0.25*linewidth_val)
                plt.plot(time, values_ens_processed[i][-3,:], color=colors[i], linestyle=linestyles[0], linewidth=0.25*linewidth_val)
                #plot the median
                plt.plot(time, values_ens_processed[i][3,:], color=colors[i], linestyle=linestyles[0], linewidth=linewidth_val)
                lines.append(plt.Line2D((0,1),(0,0),color=colors[i], linestyle=linestyles[0], linewidth=linewidth_val))
            if skill_scores and obs_values is not None and obs_time is not None:
                #print obs_time, time
                #interp_values = np.interp(obs_time, time, values[i])
                #print interp_values, obs_z, z
                ss.append(calc_skill_score(values_ens_processed[i][3,:], obs_values, 1.0, 4.0))
                #print ss
                #legend_labels[i] = legend_labels[i] + '(' + str(ss) + ')'
    if obs_values is not None and obs_time is not None:
        plt.plot(obs_time, obs_values, color='black', linewidth=linewidth_val)
        lines.append(plt.Line2D((0,1),(0,0),color='black', linewidth=linewidth_val))
        if obs_label is not None:
            legend_labels.append(obs_label)
        else:
            legend_labels.append('obs')

    #plt.legend(lines, labels, loc='best', fontsize='small')

    plt.xlabel(x_label)
    plt.ylabel(y_label)

    if multi_legend:
        first_legend = plt.legend(lines2, legend2_labels, loc='upper left', fontsize='small', bbox_to_anchor=(1.01, 1), borderaxespad=0.)
        ax = plt.gca().add_artist(first_legend)
        plt.legend(lines1, legend1_labels, loc='lower left', fontsize='small', bbox_to_anchor=(1.01, 0), borderaxespad=0.)
    else:
        if skill_scores and obs_values is not None and obs_time is not None:
            updated_legend_labels = []
            for i in range(len(values)):
                updated_legend_labels.append(legend_labels[i] + ' (' + '{0:.3f}'.format(ss[i]) + ')')
            first_legend = plt.legend(lines, updated_legend_labels, loc='best', fontsize='small')
        else:
            first_legend = plt.legend(lines, legend_labels, loc='best', fontsize='small')

    plt.savefig(filename, bbox_inches = 'tight', pad_inches = 0.1)
    plt.close()

def plot_scatter_multi(x, y, x_label, y_label, filename, x_lim=[], y_lim=[], color_index=None):
    plt.rc('text', usetex=latex_labels)

    colors = ['#e41a1c','#4daf4a','#377eb8','#984ea3','#ff7f00','#a65628','#f781bf','#ffff33']
    if color_index is None:
        color_index = 0
    else:
        if color_index > len(colors):
            raise ValueError("Invalid color index. Expected color index < " + string(len(colors)))
        else:
            #move color_index color to the front of the colors list (leave rest of list in same order)
            colors.insert(0, colors.pop(color_index))

    fig = plt.figure()

    for i in range(len(x)):
        plt.scatter(x[i], y[i], color=colors[i])
    plt.xlabel(x_label)
    plt.ylabel(y_label)

    if x_lim:
        plt.gca().set_xlim([x_lim[0],x_lim[1]])
    if y_lim:
        plt.gca().set_ylim([y_lim[0],y_lim[1]])

    plt.savefig(filename, bbox_inches = 'tight', pad_inches = 0.1)
    plt.close()
