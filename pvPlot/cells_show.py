import pyvista as pv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import colorcet as cc
import argparse
parser = argparse.ArgumentParser(prog='showMesh')

# I never could be bothered to learn python...

parser.add_argument("a", nargs='?', default="check_string_for_empty",
                    help='zone_name')
args = parser.parse_args()
if args.a == 'check_string_for_empty':
    print('I can tell that no argument was given and I can deal with that here.')
else:
    zone_name = args.a
pv.set_plot_theme("document")
voronoi_file = 'Output/cells_' + zone_name + '.vtk'
v = pv.read(voronoi_file)
label = 'Cell Diagram - Size {}'.format(v.n_points - 2)

write_file = False 
plotter = pv.Plotter(off_screen=write_file)
colour_map = cc.glasbey_light
label = 'Voronoi Cells - Size {}'.format(v.n_cells)

_ = plotter.add_mesh(v, show_edges=True, line_width=1, show_scalar_bar=False, cmap=colour_map)
_ = plotter.show_bounds(xlabel=label, ylabel='')

if write_file:
  _ = plotter.show(cpos="xy",  screenshot='Plot/' + zone_name + '.png')
else:
  _ = plotter.show(cpos="xy")


#
plotter.close()
