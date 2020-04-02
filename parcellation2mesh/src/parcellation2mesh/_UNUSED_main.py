# -*- coding: utf-8 -*-
"""
This is a skeleton file that can serve as a starting point for a Python
console script. To run this script uncomment the following lines in the
[options.entry_points] section in setup.cfg:

    console_scripts =
         fibonacci = parcellation2mesh.skeleton:run

Then run `python setup.py install` which will install the command `fibonacci`
inside your current environment.
Besides console scripts, the header (i.e. until _logger...) of this file can
also be used as template for Python modules.

Note: This skeleton file can be safely removed if not needed!
"""

import argparse
import sys
import json
import os
import nrrd
import numpy as np
import mcubes
import scipy.ndimage
import scipy.misc
from parcellation2mesh.JSONread import *
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
    jsontextfile = open(hierarchy, "r")
    jsoncontent = json.loads(jsontextfile.read())

    flat_tree = TreeIndexer.flattenTree(jsoncontent['msg'][0])
    print(json.dumps(flat_tree, indent=2))

    exit()

    search_children(jsoncontent['msg'][0])

    # reading the parcellation volume
    REFERENCE_ann, h = nrrd.read(parcellation_volume)

    # From now on, it's mainly legacy mess from the original Atlas Pipeline
    uniqueIDS = np.unique(REFERENCE_ann)
    averageZ_dict = {}
    averageY_dict = {}
    averageX_dict = {}
    uSums = []

    force_all_voxels_to_be_visible = False
    # force_all_voxels_to_be_visible = True

    continue_ = True



    #~ for uIDS in uniqueIDS:
    for ijnMain,jnMain in enumerate(region_dictionary_to_id_ALLNAME.keys()):
        # if region_dictionary_to_id_ALLNAME[jnMain] != 730:
        #     continue


        # print(region_dictionary_to_id_ALLNAME[jnMain])
        # continue

        uIDS = region_dictionary_to_id_ALLNAME[jnMain]
        if force_all_voxels_to_be_visible:
            f = os.path.join(out_dir, "mesh_"+str(uIDS)+'.obj')
        else:
            f = os.path.join(out_dir, "mesh_"+str(uIDS)+'.obj')
        if not (continue_ and os.path.isfile(f)):
            u = np.zeros(REFERENCE_ann.shape, np.float32)
            #~ u[ REFERENCE_ann==uIDS ] = 1.0
            found_at_least_one = 0

            for jn in region_dictionary_to_id_ALLNAME.keys():
                if jn.find(jnMain)>=0:
                    foundIDS = np.where(REFERENCE_ann==region_dictionary_to_id_ALLNAME[jn])
                    if foundIDS[2].shape[0]>0:
                        u[ foundIDS ] = 1.0
                        #~ voxelsXTMP += list(foundIDS[0]); voxelsYTMP += list(foundIDS[2]); voxelsZTMP += list(foundIDS[1]);
                        found_at_least_one = 1

            if found_at_least_one:
                print(region_dictionary_to_id_ALLNAME[jnMain])

                uOLD = np.copy(u)
                u = scipy.ndimage.filters.gaussian_filter( u , sigma=3.0 )

                u /= np.max(u)
                thr1 = 0.10
                u[u <thr1] = 0.0
                u[u>=thr1] = 1.0

                if force_all_voxels_to_be_visible:
                    u[uOLD>0.5] = 1.0 # forces positive voxels to be represented in the mesh

                u2 = np.zeros((u.shape[0]+2,u.shape[1],u.shape[2]), np.float32)
                u2[1:-1,:,:] = u[:,:,:]

                # Extract the 0-isosurface
                vertices, triangles = mcubes.marching_cubes(u2, 0.0)

                maxtriangles = 200000000
                if triangles.shape[0] < maxtriangles and triangles.shape[0]>3:
                    print("Exporting")
                    f = open(f, 'w')

                    for v in vertices:
                        f.write("v "+str(v[0])+" "+str(v[1])+" "+str(v[2])+" \n")

                    for t in triangles:
                        f.write("f "+str(int(t[0])+1)+" "+str(int(t[1])+1)+" "+str(int(t[2])+1)+" \n")

                    f.close()
                else:
                    print("Size too large (> "+str(maxtriangles)+" triangles) or too low...")
                    print(jnMain)
            # else:
            #     print('NOT found')

    # Smoothing the meshes, using Blender
    # command = os.path.join(folder_LIB,"Blender2.77","blender") + " -b -P mesh_smoother.py -- " + ("1 " if force_all_voxels_to_be_visible else "0 ") + folder_www
    # print(command)
    # os.system(command)







if __name__ == "__main__":
    main()
