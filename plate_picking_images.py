#! /usr/bin/env python2.7

"""
Utilities for making multi-page pdf images highlighting the colonies to pick.
 -- Weronika Patena, 2016
"""

# standard library
from __future__ import division
import sys
import string
# other packages
import matplotlib
import matplotlib.pyplot as mplt
from matplotlib.backends.backend_pdf import PdfPages
# my modules
import plotting_utilities


def plot_colony_position(axes, label_info, colony):
    mplt.axes(axes)
    xpos = int(colony[1:]) - 0.5
    ypos = 15 - string.uppercase.find(colony[0]) + 0.5
    mplt.plot(xpos, ypos, 'o', markeredgecolor='red', markersize=12, color='red')
    mplt.vlines(range(1,24), 0, 16)
    mplt.hlines(range(1,16), 0, 24)
    mplt.plot([0,0.3], [15.7,16], 'k')
    mplt.plot([0,0.5], [15.5,16], 'k')
    mplt.xlim(0,24)
    mplt.ylim(0,16)
    mplt.xticks([])
    mplt.yticks([])
    mplt.xlabel('\n' + label_info, fontsize='large')


def make_colony_picking_pdf_singles(name, position_list):
    """ Make colony-picking pdf with one image per page. 
    
    This is useful if you're displaying images on a tablet/etc.
    """
    plate_x, plate_y = 4.35, 2.93
    label_y = .4
    single_x, single_y = plate_x + 1, plate_y + label_y + 1
    with PdfPages('%s.pdf'%name) as pdf:
        for (title, well) in position_list:
            mplt.figure(figsize=(single_x, single_y))
            A = mplt.axes([.5/single_x, (.5 + label_y)/single_y, plate_x/single_x, plate_y/single_y])
            plot_colony_position(A, title, well)
            pdf.savefig()
            mplt.close()


def make_colony_picking_pdf_4perpage(name, position_list):
    """ Make colony-picking pdf with four images per page. 
    
    This is useful if you're printing the images, to save paper. 
    """
    plate_x, plate_y = 4.35, 2.93
    label_y = .4
    page_x, page_y = 11, 8.5
    padding_y = (page_y - 2*(plate_y+label_y)) / 3
    padding_x = (page_x - 2*plate_x) / 3
    xpos_1, xpos_2 = padding_x/page_x, (padding_x*2 + plate_x)/page_x
    ypos_1, ypos_2 = (padding_y+label_y)/page_y, (padding_y*2 + plate_y + label_y*2)/page_y
    def setup_page():
        mplt.figure(figsize=(page_x, page_y))
        A1 = mplt.axes([xpos_1, ypos_2, plate_x/page_x, plate_y/page_y])
        A2 = mplt.axes([xpos_2, ypos_2, plate_x/page_x, plate_y/page_y])
        A3 = mplt.axes([xpos_1, ypos_1, plate_x/page_x, plate_y/page_y])
        A4 = mplt.axes([xpos_2, ypos_1, plate_x/page_x, plate_y/page_y])
        return A1, A2, A3, A4
    position_quad_list = [position_list[x:x+4] for x in range(0, len(position_list), 4)]
    with PdfPages('%s.pdf'%name) as pdf:
        for quad in position_quad_list:
            A1, A2, A3, A4 = setup_page()
            plot_colony_position(A1, quad[0][0], quad[0][1])
            if len(quad)>1: plot_colony_position(A2, quad[1][0], quad[1][1])
            if len(quad)>2: plot_colony_position(A3, quad[2][0], quad[2][1])
            if len(quad)>3: plot_colony_position(A4, quad[3][0], quad[3][1])
            pdf.savefig()
            mplt.close()


# MAYBE-TODO make this a runnable command-line utility that takes in a mutant table, optionally sorts and checks for duplicates, and produces a sorted mutant table and a colony-picking pdf?
