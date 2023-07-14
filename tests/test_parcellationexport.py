import os, json
from rdflib import Graph, RDF
from jsonpath_ng import parse

import nrrd
import numpy as np

from blue_brain_atlas_web_exporter import parcellationexport

children = parcellationexport.children
repr_key = parcellationexport.represented
regionVolume = parcellationexport.regionVolume
regionVolumeRatio = parcellationexport.regionVolumeRatio

test_folder = os.environ["TEST_FOLDER"]
data_folder = os.path.join(test_folder, "data")
input_hierarchy = "hierarchy_l23split_first_generation.json"
input_annotation = "annotation_first_generation.nrrd"
input_annotation_path = os.path.join(data_folder, input_annotation)
output_metadata = "metadata.json"
output_hierarchy = input_hierarchy.replace(".json", "_volume.json")
output_hierarchy_path = os.path.join(data_folder, output_hierarchy)

def test_parcellationexport():
    input_hierarchy_path = os.path.join(data_folder, input_hierarchy)
    output_metadata_path = os.path.join(data_folder, output_metadata)

    parcellationexport.main_(input_hierarchy_path, input_annotation_path, "brain_region_mask", "", output_metadata_path, output_hierarchy_path, "")

def test_volumeFields():
    # Check that representedInAnnotation, regionVolume and regionVolumeRatio are correctly set
    output_hierarchy_file = open(output_hierarchy_path)
    output_hierarchy_json = json.load(output_hierarchy_file)
    represented_check(repr_key, bool, output_hierarchy_json, output_hierarchy_path)

    volume_sum = 0
    for child in output_hierarchy_json[children]:
        volume_sum += represented_check(repr_key, bool, child, output_hierarchy_path)
        offspring(child, children, represented_check, repr_key, bool, output_hierarchy_path)

    # Check that children volumes add up to parent volume
    nrrd_data, nrrd_header = nrrd.read(input_annotation_path)
    root_volume = np.count_nonzero(nrrd_data == output_hierarchy_json['id']) * parcellationexport.voxel_world_volume(nrrd_header["space directions"])

    assert (output_hierarchy_json[regionVolume]  ==  root_volume + volume_sum)


def offspring(dictionary, relation, test, test_key, test_key_type, filepath):
    for child in dictionary[relation]:
        test(test_key, test_key_type, child, filepath)
        offspring(child, relation, test, test_key, test_key_type, filepath)

def represented_check(key, key_type, dictionary, filepath):
    repres = _key_exist(key, dictionary, filepath)
    _key_type(key, key_type, dictionary, filepath)
    volume = 0
    if repres:
        volume = _key_exist(regionVolume, dictionary, filepath)
        _key_exist(regionVolumeRatio, dictionary, filepath)
    return volume

def _key_exist(key, dictionary, filepath):
    if key not in dictionary:
        although = ""
        if key != repr_key:
            although = f" although '{repr_key}' is true"
        print(f"'{key}' not found for region id {dictionary['id']} in {filepath}{although}")
        exit(1)
    else:
        return dictionary[key]

def _key_type(key, key_type, dictionary, filepath):
    if type(dictionary[key]) is not key_type:
        print(f"'{key}' is not of type {key_type} for region id {dictionary['id']} in {filepath}")
        print(f"'{key}': {dictionary[key]}")
        exit(1)

