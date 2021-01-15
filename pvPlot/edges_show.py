import pyvista as pv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import colorcet as cc
import argparse
parser = argparse.ArgumentParser(prog='showMesh')


parser.add_argument("a", nargs='?', default="check_string_for_empty",
                    help='zone_name')
args = parser.parse_args()
if args.a == 'check_string_for_empty':
    print('I can tell that no argument was given and I can deal with that here.')
else:
    zone_name = args.a

pv.set_plot_theme("document")
voronoi_file = 'Output/edges_' + zone_name + '.vtk'
v = pv.read(voronoi_file)

label = 'Triangulation & Voronoi Edges - Size {}'.format(v.n_cells - 2)
#v["Labels"] = ["{}".format(v.active_scalars[i]) for i in range(v.n_points)]

write_file = False
plotter = pv.Plotter(off_screen=write_file)
colour_map = cc.glasbey_dark

#_ = plotter.add_point_labels(v, "Labels", point_size=5, font_size=10)
_ = plotter.add_mesh(v, show_scalar_bar=False, cmap=colour_map)
#_ = plotter.add_mesh(v, show_scalar_bar=False, cmap=colour_map, color='yellow')
#_ = plotter.show_bounds(xlabel=label, ylabel='')

if write_file:
  _ = plotter.show(cpos="xy", screenshot='Plot/' + zone_name + '.png')
else:
  _ = plotter.show(cpos="xy")


#
plotter.close()
