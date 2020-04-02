# -*- coding: utf-8 -*-
"""

"""

import argparse
import sys
import json
import os
import nrrd
import numpy as np
# import mcubes
from skimage import measure
import optimesh


from parcellation2mesh import __version__

__author__ = "jonathanlurie"
__copyright__ = "jonathanlurie"
__license__ = "mit"



def parse_args(args):
    """Parse command line parameters

    Args:
      args ([str]): command line parameters as list of strings

    Returns:
      :obj:`argparse.Namespace`: command line parameters namespace
    """
    parser = argparse.ArgumentParser(
        description="Just a Fibonacci demonstration")
    parser.add_argument(
        "--version",
        action="version",
        version="parcellation2mesh {ver}".format(ver=__version__))

    parser.add_argument(
        "--hierarchy",
        dest="hierarchy",
        required=True,
        metavar="<FILE PATH>",
        help="The hierarchy JSON file, sometimes called 1.json")

    parser.add_argument(
        "--parcellation-volume",
        dest="parcellation_volume",
        required=True,
        metavar="<FILE PATH>",
        help="The NRRD parcellation volume file")

    parser.add_argument(
        "--out-dir",
        dest="out_dir",
        required=True,
        metavar="<DIRECTORY PATH>",
        help="The output directory to create the OBJ mesh files")


    return parser.parse_args(args)


def exportObj(vertices, triangles, filepath):
    f = open(filepath, 'w')

    for v in vertices:
        f.write("v "+str(v[0])+" "+str(v[1])+" "+str(v[2])+" \n")

    for t in triangles:
        # f.write("f "+str(int(t[0])+1)+" "+str(int(t[1])+1)+" "+str(int(t[2])+1)+" \n")
        f.write("f "+str(int(t[2])+1)+" "+str(int(t[1])+1)+" "+str(int(t[0])+1)+" \n")

    f.close()


def main():
    """Main entry point allowing external calls

    Args:
      args ([str]): command line parameter list
    """
    args = parse_args(sys.argv[1:])

    hierarchy = args.hierarchy
    parcellation_volume = args.parcellation_volume
    out_dir = args.out_dir

    # create out_dir if inexistant
    try:
        os.makedirs(out_dir)
    except FileExistsError as e:
        pass

    # loading json annotation


    # reading the parcellation volume
    print("Reading NRRD...")
    parcellation_data, parcellation_header = nrrd.read(parcellation_volume)
    # present_regions = np.unique(parcellation_data).tolist()

    print("Prepare outer hull volume...")
    outer_hull = np.zeros_like(parcellation_data, dtype = "uint8")

    outer_hull[parcellation_data != 0] = 255

    print("Marching cube...")
    # vertices, triangles = mcubes.marching_cubes(outer_hull, 0)
    # vertices, faces, normals, values = measure.marching_cubes(outer_hull, 0)
    # vertices, faces, normals, values = measure.marching_cubes_lewiner(p, threshold, step_size=step_size, allow_degenerate=True)
    # vertices, triangles = measure.marching_cubes_classic(outer_hull)
    vertices, triangles, normals, values = measure.marching_cubes_lewiner(outer_hull)

    print("Smoothing...")
    # vertices_smooth, triangles_smooth = optimesh.cpt.fixed_point_uniform(
    vertices_smooth, triangles_smooth = optimesh.odt.fixed_point_uniform(
    # vertices_smooth, triangles_smooth = optimesh.cvt.quasi_newton_uniform_full(
        vertices, triangles, 1.0e-2, 100, verbose=True
        # implicit_surface=Sphere(),
        # step_filename_format="out{:03d}.vtk"
    )

    print("Export OBJ...")
    rough_mesh_filepath = os.path.join(out_dir, "outer_hull_smooth_rough_lewiner_ccfv3_optimesh.obj")
    exportObj(vertices_smooth, triangles_smooth, rough_mesh_filepath)


if __name__ == "__main__":
    main()
