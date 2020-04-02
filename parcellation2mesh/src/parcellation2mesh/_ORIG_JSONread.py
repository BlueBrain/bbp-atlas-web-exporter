"""tools to deal with brain hierarchy from the AIBS"""
import numpy as np

id_to_region_dictionary = {}
id_to_region_dictionary_ALLNAME = {}  # id to complete name
region_dictionary_to_id = {}  # name to id
region_dictionary_to_id_ALLNAME = {}  # complete name to id

region_dictionary_to_id_ALLNAME_parent = {
}  # complete name to complete name parent
region_dictionary_to_id_parent = {}  # name to complete name parent
allname2name = {}  # complete name to name
name2allname = {}  # name to complete name
region_keys = []  # list of regions names
regions_ALLNAME_list = []  # list of complete regions names
is_leaf = {}  # full name to int (! if is leaf, else 0)
id_to_color = {}  # region id to color in RGB
region_to_color = {}  # complete name to color in RGB


def return_ids_containing_str_list(str_list):
    """
    Retrieve the list of region id which complete name contains all the keywords in str_list.

    Arguments:
        str_list: List of keyword that the region complete name.
    Returns:
        List of region id matching condition
    """

    id_list = []
    for kk in id_to_region_dictionary_ALLNAME:
        region_is_in = True
        for str1 in str_list:
            if (id_to_region_dictionary_ALLNAME[kk].lower()).find(str1.lower(
            )) < 0:  # if any of the regions is not there, do not take
                region_is_in = False
                break
        if region_is_in:
            id_list.append(kk)
    return id_list


def hex_to_rgb(value):
    """
    Converts a Hexadecimal color into its RGB value counterpart.

    Arguments:
        value: string hexadecimal color to convert.
    Returns:
        List of the Red, Green, and Blue components of the color
    """

    value = value.lstrip('#')
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))


def search_children(object_,
                    lastname_ALL="",
                    lastname="",
                    darken=True):
    """
    Explores the hierarchy dictionary to extract its brain regions and fills external dictionaries.
    Arguments
        object_: dictionary of regions properties. See
        https://bbpteam.epfl.ch/documentation/projects/voxcell/latest/atlas.html#brain-region-hierarchy
        lastname_ALL: complete name of the parent of the current brain region
        lastname: name of the parent of the current brain region
        darken: if True, darkens the region colors too high
    """

    regions_ALLNAME_list.append(lastname_ALL + "|" + object_["name"])
    name2allname[object_["name"]] = lastname_ALL + "|" + object_["name"]
    allname2name[lastname_ALL + "|" + object_["name"]] = object_["name"]
    id_to_region_dictionary[object_["id"]] = object_["name"]
    id_to_region_dictionary_ALLNAME[
        object_["id"]] = lastname_ALL + "|" + object_["name"]
    region_dictionary_to_id[object_["name"]] = object_["id"]
    region_dictionary_to_id_ALLNAME[lastname_ALL + "|" +
                                    object_["name"]] = object_["id"]
    region_dictionary_to_id_ALLNAME_parent[lastname_ALL + "|" +
                                           object_["name"]] = lastname_ALL
    region_dictionary_to_id_parent[object_["name"]] = lastname
    clrTMP = np.float32(
        np.array(list(hex_to_rgb(object_["color_hex_triplet"]))))
    if np.sum(clrTMP) > 255.0 * 3.0 * 0.75 and darken:
        clrTMP *= 255.0 * 3.0 * 0.75 / np.sum(clrTMP)
    region_to_color[lastname_ALL + "|" + object_["name"]] = list(clrTMP)
    id_to_color[object_["id"]] = list(clrTMP)
    region_keys.append(object_["name"])
    try:
        is_leaf[lastname_ALL + "|" + object_["name"]] = 1
        # ~ region_dictionary_to_id_ALLNAME_child[  lastname_ALL+"|"+object_["name"] ] = children
        # ~ id_children[object_["id"]] = object_["children"]
        for children in object_["children"]:
            search_children(children,
                            lastname_ALL + "|" + object_["name"],
                            object_["name"], darken=darken)
            is_leaf[lastname_ALL + "|" + object_["name"]] = 0
    except KeyError:
        print("No children of object")


dict_corrections = {}
old_regions_layer23 = [
    41, 113, 163, 180, 201, 211, 219, 241, 251, 269, 288, 296, 304, 328, 346,
    412, 427, 430, 434, 492, 556, 561, 582, 600, 643, 657, 667, 670, 694, 755,
    806, 821, 838, 854, 888, 905, 943, 962, 965, 973, 1053, 1066, 1106, 1127,
    12994, 182305697
]
for reg in old_regions_layer23:
    dict_corrections[reg] = [reg + 20000, reg + 30000]

# Change of id when L2 and L2/3 existed
dict_corrections[195] = [20304]
dict_corrections[747] = [20556]
dict_corrections[524] = [20582]
dict_corrections[606] = [20430]

inv_corrections = {}
for k, v in dict_corrections.items():
    for conv in v:
        inv_corrections[conv] = k
