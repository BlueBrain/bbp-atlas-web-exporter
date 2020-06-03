# -*- coding: utf-8 -*-
"""

"""

import argparse
import sys
import json
import os
import nrrd
import numpy as np
from skimage import measure
import parcellation2mesh.TreeIndexer as TreeIndexer


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
        description="Export volumetric parcellations to meshes.")
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

    # reading the nrrd file
    nrrd_data, nrrd_header = nrrd.read(parcellation_volume)

    # loading json annotation
    jsoncontent = json.loads(open(hierarchy, "r").read())

    # sometimes, the 1.json has its content in a "msg" sub prop (the original vesion has).
    # and some other versions don't. Here we deal with both
    if "msg" in jsoncontent:
        flat_tree = TreeIndexer.flattenTree(jsoncontent['msg'][0])
    else:
        flat_tree = TreeIndexer.flattenTree(jsoncontent)

    # print(json.dumps(flat_tree, indent=2))

    total_region = len(flat_tree)
    region_counter = 0

    # For each region, we create a mask that contains all the sub regions
    for region_id in flat_tree:
        region_counter += 1
        region_node = flat_tree[region_id]

        print("{}/{} - [{}] {}".format(region_counter, total_region, region_id, flat_tree[region_id]["name"]))

        region_mask = np.zeros_like(nrrd_data, dtype = "uint8")

        # masking the current region
        print("Creating region mask...")
        region_mask[nrrd_data == region_id] = 255
        subregion_counter = 0
        total_subregions = len(region_node["_descendants"])
        for child_id in region_node["_descendants"]:
            subregion_counter += 1
            print("Grouping subregions {}/{}".format(subregion_counter, total_subregions), end="\r")
            # masking the current region
            region_mask[nrrd_data == child_id] = 255

        print()

        # if the mask is all black, then there is no mesh to build
        if not np.any(region_mask):
            continue

        # Creating the mesh with the marching cube
        print("Marching cube...")
        vertices, triangles, normals, values = measure.marching_cubes_lewiner(region_mask)

        # Exporting the mesh as OBJ file
        print("Export OBJ...")
        rough_mesh_filepath = os.path.join(out_dir, str(region_id) + ".obj")
        exportObj(vertices, triangles, rough_mesh_filepath)


if __name__ == "__main__":
    main()
