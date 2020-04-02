import sys
import os
from os.path import *
import numpy as np
import nrrd
import json
import scipy.ndimage
import scipy.misc
from JSONread import *
import mcubes

if len(sys.argv)<2:
    print('The data folder path is required')
    sys.exit(1)
DATA_FOLDER = sys.argv[-1]
VOLUME_RESOLUTION = 25.0
folder_LIB = "../../LIB/"

folder_www = join(DATA_FOLDER, "website_export")
# loading json annotation
jsontextfile = open(join(DATA_FOLDER, "allen_data", "1fixed.json"), "r")
jsoncontent = json.loads(jsontextfile.read())
search_children(jsoncontent['msg'][0])


REFERENCE_ann, h = nrrd.read( join(DATA_FOLDER, "allen_data", "annotation_" + str(int(VOLUME_RESOLUTION)) + "fixed.nrrd") )

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
    print(region_dictionary_to_id_ALLNAME[jnMain])
    uIDS = region_dictionary_to_id_ALLNAME[jnMain]
    if force_all_voxels_to_be_visible:
        folder_mesh = join(folder_www, "marching_cubes", "meshesT")
        os.system("mkdir -p "+ folder_mesh)
        f = join(folder_mesh, "mesh_"+str(uIDS)+'.obj')
    else:
        folder_mesh = join(folder_www, "marching_cubes", "meshesMS")
        os.system("mkdir -p "+folder_mesh)
        f = join(folder_mesh, "mesh_"+str(uIDS)+'.obj')
    if not (continue_ and isfile(f)):
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

command = join(folder_LIB,"Blender2.77","blender") + " -b -P mesh_smoother.py -- " + ("1 " if force_all_voxels_to_be_visible else "0 ") + folder_www
print(command)
os.system(command)
