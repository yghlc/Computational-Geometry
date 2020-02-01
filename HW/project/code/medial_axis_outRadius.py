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
import matplotlib.pyplot as plt
import os

# for the case running in other folder
sys.path.insert(0,os.path.expanduser('~/codes/PycharmProjects/yghlc_Computational-Geometry/HW/lib'))
# local modules
sys.path.insert(0, '../../lib/')
from polygon_medial_axis import compute_polygon_medial_axis,\
        plot_polygon_medial_axis


def main():
    if len(sys.argv) < 2:
        print('Usage: >> python {} <polygon_file>'.format(sys.argv[0]))
        sys.exit(1)
    polygon = np.loadtxt(sys.argv[1])
    h = 0.1
    if len(sys.argv)==3:
        h = float(sys.argv[2])
    fig, ax = plt.subplots(figsize=(8, 8))
    medial_axis, radiuses = compute_polygon_medial_axis(polygon, h=h)


    # medial_axis = compute_polygon_medial_axis(polygon, h=0.1)
    # plot_polygon_medial_axis(polygon, medial_axis, ax=ax)
    # ax.axis('equal')
    # ax.set_title('Medial Axis')
    # plt.show()
    # plt.savefig('fig.jpg')

    # save radius to file
    with open('save_medial_axis_radius.txt','w') as f_obj:
        for ((x1,y1),(x2,y2)),(r1, r2) in zip(medial_axis, radiuses):
            f_obj.writelines('%f  %f  %f %f  %f  %f\n'%(x1, y1, x2, y2,r1, r2))


if __name__ == "__main__":
    main()
