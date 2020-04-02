import bpy
import os
from os.path import *
import sys

argv = sys.argv
force_all_voxels_to_be_visible = argv[argv.index("--") + 1] == "1"
data_folder = argv[argv.index("--") + 2]
if force_all_voxels_to_be_visible:
    meshes_folder = 'meshesT'
else:
    meshes_folder = 'meshesMS'

folder_ = join(data_folder, "marching_cubes", meshes_folder)
output_folder = join(data_folder, "COMDATA", meshes_folder)
os.system("mkdir -p "+output_folder)
print(folder_)
filelist = os.listdir(folder_)
for ff in filelist:
    # print(ff)
    if ff[:6]!="smooth" and ff[:8]!="decimate":
        meshIDstr = (ff.split("_")[-1]).split(".")[0]
        bpy.ops.import_scene.obj(filepath=join(folder_,ff), filter_glob="*.obj;*.mtl", use_edges=True, use_smooth_groups=True, use_split_objects=True, use_split_groups=True, use_groups_as_vgroups=False, use_image_search=True, split_mode='ON', global_clamp_size=0, axis_forward='-Z', axis_up='Y')
        mod = bpy.data.objects["mesh_"+meshIDstr].modifiers.new(name='smooth', type='SMOOTH')
        mod.factor = 1.0
        mod.iterations = 10
        if len(bpy.data.objects["mesh_"+meshIDstr].data.vertices)>400:
            mod2 = bpy.data.objects["mesh_"+meshIDstr].modifiers.new(name='decimate', type='DECIMATE')
            mod2.decimate_type = 'COLLAPSE'
            mod2.ratio = 0.1
        bpy.context.scene.objects.active = bpy.data.objects["mesh_"+meshIDstr]
        bpy.ops.export_scene.obj(filepath=output_folder+"/decimated_smoothed_"+ff, use_selection=True)
        bpy.ops.object.delete()

bpy.ops.wm.quit_blender()
