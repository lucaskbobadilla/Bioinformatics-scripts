#!/bin/python


# Modules necessary for running this script

import cairo
import math
import re 
import numpy as np
import sys
import argparse


#
# The img dictionary will hold image-specific dimensions. 
#   

img = {}
img['height'] = 800
img['width'] = 800
img['center_x'] = img['width'] / 2.0 
img['center_y'] = img['height'] / 2.0 
img['font_size'] = 16
img['radius_rna_1']   = img['width'] * 0.27
img['radius_rna_2']   = img['width'] * 0.30
img['radius_fst_1']   = img['width'] * 0.31
img['radius_fst_2']   = img['width'] * 0.34 
img['radius_group_label'] = img['width'] * 0.37
img['tick_outer_radius'] = img['width'] * 0.35
img['tick_inner_radius'] = img['width'] * 0.25
img['tick_label_radius'] = img['width'] * 0.235
img['fst_legend_x'] = img['width'] * 0.93
img['rna_legend_x'] = img['width'] * 0.87
img['legend_width'] = 30
img['legend_start_y'] = img['height'] * 0.05
img['legend_height'] = img['legend_width'] * 7

# get values from dict 
center_x = img['center_x']
center_y = img['center_y']
fontsize = img['font_size']
radius_group_label = img['radius_group_label']
radius_1 = img['radius_rna_1']
radius_2 = img['radius_rna_2']
radius_3 = img['radius_fst_1']
radius_4 = img['radius_fst_2']
tick_outer_radius = img['tick_outer_radius']
tick_inner_radius = img['tick_inner_radius']
fst_legend_x, fst_legend_y = img['fst_legend_x'], img['legend_start_y']
rna_legend_x, rna_legend_y = img['rna_legend_x'], img['legend_start_y']
legend_width = img['legend_width']
legend_height = img['legend_height']
tick_label_radius = img['tick_label_radius']

#set figure start position
start = 2.0
end = 0

#
# get data sets
#

p = argparse.ArgumentParser() #parse arguments

# name arguments parsing
p.add_argument('-e', type=str, required=True, help = "This argument is the expression file path")
p.add_argument('-f', type=str, required=True, help = "This argument is the fst stats file path")

args = p.parse_args() # get arguments

# open files
expression_df = open(args.e, 'r')
fst_df = open(args.f, 'r')

#
# get fst values into dict
#
fst_stats = {}
for line in fst_df:
    line = line.strip("\n")
    if line[0] != "#":
        line_break = line.split("\t")
        scaf_match = re.match("^scaffold",line_break[4])
        if scaf_match:
            continue
        elif 'bp' and 'stat' not in fst_stats.get(line_break[4],{}):
            fst_stats.setdefault(line_break[4], {})['bp'] = []
            fst_stats.setdefault(line_break[4], {})['stat'] = []
        else:
            fst_stats[line_break[4]]['bp'].append(int(line_break[5]))
            fst_stats[line_break[4]]['stat'].append(float(line_break[10]))

#
# get expression values into dict 
#
rna_stats = {}
for line in expression_df:
    line = line.strip("\n")
    if line[0:3] != "Gene":
        line_break = line.split("\t")
        group_match = re.match("^group",line_break[1])
        if not group_match:
            continue
        elif 'bp' and 'stat' not in rna_stats.get(line_break[1],{}):
            rna_stats.setdefault(line_break[1], {})['bp'] = []
            rna_stats.setdefault(line_break[1], {})['stat'] = []
        elif line_break[5] == "TRUE":
            rna_stats[line_break[1]]['bp'].append(int(line_break[2]))
            rna_stats[line_break[1]]['stat'].append(float(line_break[6]))


#
# Convert a radius and a span of degrees into X, Y coordinates 
#
def get_x_y_coordinates(center_x, center_y, degree, radius):
    if degree <= 90:
        theta = float(degree)
        opp_side = radius * math.sin(math.radians(theta))
        adj_side = radius * math.cos(math.radians(theta))
        x = center_x + adj_side
        y = center_y + opp_side
    elif degree <= 180:
        theta = float(degree - 90.0)
        opp_side = radius * math.sin(math.radians(theta))
        adj_side = radius * math.cos(math.radians(theta))
        x = center_x - opp_side
        y = center_y + adj_side
    elif degree <= 270:
        theta = float(degree - 180.0)
        opp_side = radius * math.sin(math.radians(theta))
        adj_side = radius * math.cos(math.radians(theta))
        x = center_x - adj_side
        y = center_y - opp_side
    else:
        theta = float(degree - 270.0)
        opp_side = radius * math.sin(math.radians(theta))
        adj_side = radius * math.cos(math.radians(theta))
        x = center_x + opp_side
        y = center_y - adj_side
    return (x, y)

#
# Chromosome lengths
#
chrm_len = {'I' : 28185914, 'II' : 23295652, 'III' : 16798506, 'IV' : 32632948,
'IX' : 20249479, 'V' : 12251397, 'VI' : 17083675, 'VII' : 27937443,
'VIII' : 19368704, 'X' : 15657440, 'XI' : 16706052, 'XII' : 18401067,
'XIII' : 20083130, 'XIV' : 15246461, 'XIX' : 20240660, 'XV' : 16198764,
'XVI' : 18115788, 'XVII' : 14603141, 'XVIII' : 16282716, 'XX' : 19732071,
'XXI' : 11717487}

#
# define total genome size
#
total_genome = sum(chrm_len.values())

#
#  define function for chromosome scaling
#
def chrom_scale(chrm_length,total_genome, groups):
    margin_groups = 2.0*groups # margin of 2 to fit all
    chrom_size = round((360.0-margin_groups) * (float(chrm_length)/float(total_genome)),10)
    return chrom_size

#
# define function to get colour
#
def  pick_colour(stat, name):
    if name == "fst":
        return (abs(stat), 0,abs(stat)) # transform all stat value int absolute values to take any negative value out
    elif name == "rna":
        # convert to a log scale and multiple it by 0.1 to fit into a 0 to 1 scale and use absolute values only
        stat = abs(stat*0.01) 
        return (0.5+stat, 0.8-stat,0)
    else:
        return "name need to be fst or rna"

#
# define function for chromosome drawing 
#
def chrm_draw(start, end, inner_radius, outer_radius,r,g,b):
    #first x and y
    x1, y1 = get_x_y_coordinates(center_x, center_y, start, inner_radius)
    ctx.move_to(x1,y1)
    ctx.arc(center_x, center_y,inner_radius,math.radians(start), math.radians(end))
    x2, y2 = get_x_y_coordinates(center_x, center_y, end, outer_radius)
    ctx.line_to(x2,y2)
    ctx.arc_negative(center_x, center_y, outer_radius, math.radians(end), math.radians(start))
    ctx.close_path()
    ctx.set_source_rgb(r,g,b)
    ctx.fill_preserve()
    ctx.set_source_rgb(105/255, 105/255, 105/255)
    ctx.set_line_width(2)
    ctx.stroke()

#
# define function to draw coordinates 
#
def coord_line(group_length, start, tick_length,total_genome):
    step_tick_radian = chrom_scale(tick_length, total_genome, 21)
    tick_loc = np.arange(0,group_length,int(tick_length))
    for i in tick_loc:
        tick_position = i / tick_length
        tick_radian = start + (step_tick_radian * tick_position)
        if tick_position % 2 == 0:
            xIn, yIn = get_x_y_coordinates(center_x, center_y, tick_radian, tick_inner_radius)
            ctx.move_to(xIn, yIn)
            xOut, yOut = get_x_y_coordinates(center_x, center_y, tick_radian, tick_outer_radius)
            ctx.line_to(xOut, yOut)
            ctx.set_source_rgb(0,0,0)
            ctx.set_line_width(1)
            ctx.stroke()
            if tick_position != 0:
                name = str(round(5 * tick_position)) + 'Mb'
                label_tick(name, tick_radian)
        else:
            xIn, yIn = get_x_y_coordinates(center_x, center_y, tick_radian, tick_inner_radius+5)
            ctx.move_to(xIn, yIn)
            xOut, yOut = get_x_y_coordinates(center_x, center_y, tick_radian, tick_outer_radius-3)
            ctx.line_to(xOut, yOut)
            ctx.set_source_rgb(0,0,0)
            ctx.set_line_width(.5)
            ctx.stroke()


#
# define chrmom naming label function
#
def group_label(name, radian_begin, radian_end):
    ctx.set_font_size(fontsize)
    ctx.set_source_rgb(0,0,0)
    text = name
    extents = ctx.text_extents(text)
    text_width = extents[2]
    text_height = extents[3]
    radian_group_label = (radian_begin + radian_end)/2
    group_x, group_y = get_x_y_coordinates(center_x - (center_x * 0.020), center_y + (center_y * 0.015), radian_group_label,radius_group_label)
    ctx.move_to(group_x, group_y)
    ctx.show_text(text)

#
# define stat_line function for adding statistics values into the arcs
#
def stats_line(name,value, bp, inner_radius, outer_radius):
    bp_angle = chrom_scale(start+bp,total_genome,21)
    bp_position = start + bp_angle
    x1bp, y1bp = get_x_y_coordinates(center_x, center_y, bp_position, inner_radius+2)
    ctx.move_to(x1bp, y1bp)
    x2bp, y2bp = get_x_y_coordinates(center_x, center_y, bp_position, outer_radius-2)
    ctx.line_to(x2bp,y2bp)
    colours = pick_colour(value,name)
    ctx.set_source_rgb(colours[0], colours[1], colours[2])
    ctx.set_line_width(.5)
    ctx.stroke()

#
# Draw legend box function
#
def draw_legend_box(variable, x, y, legend_width, legend_height):
    #get rectangle coord
    ctx.move_to(x, y)
    ctx.line_to(x, y + legend_height+1)
    ctx.line_to(x + legend_width, y + legend_height+1)
    ctx.line_to(x + legend_width, y)
    ctx.close_path()
    ctx.set_line_width(1)
    if variable == "fst":
        ctx.set_source_rgb(0,0,0)
        ctx.set_line_width(1.2)
    elif variable == "rna":
        ctx.set_source_rgb(0,0,0)
        ctx.set_line_width(1.2)
    ctx.fill()

#
# add legend label function
#
def legend_label(name, legend_start_x, legend_start_y):
    ctx.set_font_size(fontsize - 6)
    ctx.set_source_rgb(0,0,0)
    text = name
    extents = ctx.text_extents(text)
    text_width = extents[2]
    text_height = extents[3]
    ctx.move_to(legend_start_x + legend_width/2 - 6, legend_start_y - 5)
    ctx.show_text(text)

#
# add legend numbers function
#
def legend_number(value_type):
    if value_type == 'fst':
        levels = ['1.0', '0.5', '0.0']
        legend_x = fst_legend_x 
        legend_y = fst_legend_y
        for item in levels:
            text = item
            ctx.set_font_size(fontsize - 6)
            ctx.move_to(legend_x-19, legend_y)
            ctx.show_text(text)
            legend_y = legend_y + (legend_height/2)
    if value_type == 'rna':
        levels = ['1.0', '0.5', '0.0']
        legend_x = rna_legend_x 
        legend_y = rna_legend_y
        for item in levels:
            text = item
            ctx.set_font_size(fontsize - 6)
            ctx.move_to(legend_x-19, legend_y)
            ctx.show_text(text)
            legend_y = legend_y + (legend_height/2)

#
# define legend levels function
#
def legend_key(variable):
    if variable == "fst":
        fst_values_leg = np.arange(1.0001,0.0,-0.0001)
        x = fst_legend_x
        y = fst_legend_y
        for value in fst_values_leg:
            ctx.move_to(x, y)
            ctx.line_to(x+legend_width-6, y)
            colours = pick_colour(value,'fst')
            ctx.set_source_rgb(colours[0], colours[1], colours[2])
            ctx.set_line_width(.5)
            ctx.stroke()
            increment = legend_height/len(fst_values_leg)
            y = y + increment
    elif variable == "rna":
        rna_values_leg = np.arange(100,0.0,-0.01)
        x = rna_legend_x
        y = rna_legend_y
        for value in rna_values_leg:
            ctx.move_to(x, y)
            ctx.line_to(x+legend_width-6, y)
            colours = pick_colour(value,'rna')
            ctx.set_source_rgb(colours[0], colours[1], colours[2])
            ctx.set_line_width(.5)
            ctx.stroke()
            increment = legend_height/len(rna_values_leg)
            y = y + increment

#
# define tick label function
#
def label_tick(number, radian):
    text = number
    extents = ctx.text_extents(text)
    text_width = extents[2]
    text_height = extents[3]
    radian_label = radian - 1
    x, y = get_x_y_coordinates(center_x - (center_x * 0.025), center_y + (center_y * 0.005), radian_label, tick_label_radius)
    ctx.move_to(x, y)
    ctx.set_source_rgb(0, 0, 0)
    ctx.set_font_size(img['font_size']/2)
    ctx.show_text(text)

######### Draw figure ###########

with cairo.PDFSurface("Circle_plot.pdf", img['height'], img['width']) as surface:
    ctx = cairo.Context(surface)
    for group, length in chrm_len.items():
        angle = chrom_scale(length, total_genome, 21)
        end = angle + start

        # draw chromosome arcs
        coord_line(length, start, 5000000, total_genome)
        chrm_draw(start, end, radius_3, radius_4, 0,0,0)
        chrm_draw(start,end, radius_1, radius_2, 0.5,0.8,0)

        #draw legend boxes
        draw_legend_box('fst',fst_legend_x-2, fst_legend_y-2, legend_width-2, legend_height+3)
        legend_key('fst')
        draw_legend_box('rna',rna_legend_x-2, rna_legend_y-2, legend_width-2, legend_height+3)
        legend_key('rna')

        #get arrays for stats and basepairs of fst
        group_fst_stat = fst_stats["group" + group]['stat']
        group_fst_bp = fst_stats["group" + group]['bp']
        for i in range(len(group_fst_bp)):
            value = group_fst_stat[i]
            bp = group_fst_bp[i]
            stats_line('fst',value,bp,radius_3, radius_4)

        #get arrays for fst stats and basepairs of expression
        group_rna_stat = rna_stats["group" + group]['stat']
        group_rna_bp = rna_stats["group" + group]['bp']

        # loop over arrays to add statistics into arcs
        for i in range(len(group_rna_bp)):
            value = group_rna_stat[i]
            bp = group_rna_bp[i]
            stats_line('rna',value,bp,radius_1, radius_2)
        
        # add group and legends
        group_label(group, start, end)
        start = end + 2
        legend_label('Fst',fst_legend_x-5, fst_legend_y)
        legend_number('fst')
        legend_label('FC',rna_legend_x-4, rna_legend_y)
        legend_number('rna')

#
#close files
#
expression_df.close()
fst_df.close()
