#!/usr/bin/env python3
# -*- coding: utf-8 -*-
##################################
# University of Wisconsin-Madison
# Author: Yaqi Zhang
##################################
"""
compute medial axis of a simple polygon
"""

# standard library
import sys

# third party library
import numpy as np

import os

# for the case running in other folder
sys.path.insert(0,os.path.expanduser('~/codes/PycharmProjects/yghlc_Computational-Geometry/HW/lib'))
# local modules
sys.path.insert(0, '../../lib/')
from polygon_medial_axis import compute_polygon_medial_axis,\
        plot_polygon_medial_axis

import re

def main():
    if len(sys.argv) < 2:
        print('Usage: >> python {} <polygon_file>'.format(sys.argv[0]))
        sys.exit(1)
    polygon_txt = sys.argv[1]
    polygon = np.loadtxt(polygon_txt)
    h = 0.1
    if len(sys.argv)==3:
        h = float(sys.argv[2])

    medial_axis, radiuses = compute_polygon_medial_axis(polygon, h=h)

    # avoid import matplotlib if don't need it, or it will ask for graphic environment, and make window loses focus
    # import matplotlib.pyplot as plt
    # fig, ax = plt.subplots(figsize=(8, 8))
    # medial_axis = compute_polygon_medial_axis(polygon, h=0.1)
    # plot_polygon_medial_axis(polygon, medial_axis, ax=ax)
    # ax.axis('equal')
    # ax.set_title('Medial Axis')
    # plt.show()
    # plt.savefig('fig.jpg')

    save_txt = 'save_medial_axis_radius.txt'
    proc_id_str = re.findall(r'\d+',os.path.basename(polygon_txt))
    if len(proc_id_str) == 1:
        proc_id = int(proc_id_str[0])
        save_txt = 'save_medial_axis_radius_%d.txt'%proc_id

    # save radius to file
    with open(save_txt,'w') as f_obj:
        for ((x1,y1),(x2,y2)),(r1, r2) in zip(medial_axis, radiuses):
            f_obj.writelines('%f  %f  %f %f  %f  %f\n'%(x1, y1, x2, y2,r1, r2))


if __name__ == "__main__":
    main()
