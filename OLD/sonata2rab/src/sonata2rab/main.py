import sys
import argparse
import h5py
import numpy as np
import json
import randomaccessbuffer as rab
import pandas as pd
from scipy.spatial.transform import Rotation
import sonata2rab.TreeIndexer as TreeIndexer
from sonata2rab import __version__


GROUP_0_PATH = "nodes/atlas_cells/0"
EXPECTED_MODALITY_NAMES = {
    "position": ["x", "y", "z"],
    "orientation": ["orientation_w", "orientation_x", "orientation_y", "orientation_z"],
    "region": ["region_id"],
    "cell_type": ["cell_type"],
  }


CELL_TYPE_ID_TO_CELL_TYPE_NAME = {
    "lut": [
      "nsg:ExcitatoryNeuron",
      "nsg:InhibitoryNeuron",
      "nsg:Microglia",
      "nsg:Astrocyte",
      "nsg:Oligodendrocyte",
    ]
  }


def check_positive_int(value):
    ivalue = int(value)
    if ivalue <= 0:
        raise argparse.ArgumentTypeError("%s is an invalid positive int value" % value)
    return ivalue


def parse_args(args):
  """Parse command line parameters

  Args:
    args ([str]): command line parameters as list of strings

  Returns:
    :obj:`argparse.Namespace`: command line parameters namespace
  """
  parser = argparse.ArgumentParser(
    description="Export metadata (.json) and RandomAccessBuffer (.rab) files from SONATA (.h5) cell position files.")
  parser.add_argument(
    "--version",
    action="version",
    version="sonata2rab {ver}".format(ver=__version__))

  parser.add_argument(
    "--hierarchy",
    dest="hierarchy_filepath",
    required=True,
    metavar="<FILE PATH>",
    help="The hierarchy JSON file, sometimes called 1.json")

  parser.add_argument(
    "--sonata",
    dest="sonata_filepath",
    required=True,
    metavar="<FILE PATH>",
    help="Path to the SONATA cell position file")

  parser.add_argument(
    "--downsampling-factor",
    dest="downsampling_factor",
    required=False,
    metavar="<INTEGER>",
    default=1,
    type=check_positive_int,
    help="Factor to downsample data. Default: 1 (no downsampling). A ratio of 100 will keep only one value every 100.")

  parser.add_argument(
    "--out-postion",
    dest="out_position",
    required=True,
    metavar="<FILE PATH>",
    help="Path to write the RandomAccessBuffer file for cell positions")

  parser.add_argument(
    "--out-cell-type",
    dest="out_cell_type",
    required=True,
    metavar="<FILE PATH>",
    help="Path to write the RandomAccessBuffer file for cell types")

  parser.add_argument(
    "--out-orientation",
    dest="out_orientation",
    required=True,
    metavar="<FILE PATH>",
    help="Path to write the RandomAccessBuffer file for cell orientations")

  parser.add_argument(
    "--out-region",
    dest="out_region",
    required=True,
    metavar="<FILE PATH>",
    help="Path to write the RandomAccessBuffer file for cell regions")

  parser.add_argument(
    "--out-metadata",
    dest="out_metadata",
    required=True,
    metavar="<FILE PATH>",
    help="Path to write the JSON file for region-based metadata")


  return parser.parse_args(args)
  

# Check that all the expected datasets are present in the SONATA file
def check_sonata(sonata_instance, raise_exception = True):
  group_0 = sonata_instance[GROUP_0_PATH]
  dataset_names = list(group_0.keys())

  for modality_name in EXPECTED_MODALITY_NAMES:
    for i in range(0, len(EXPECTED_MODALITY_NAMES[modality_name])):
      name = EXPECTED_MODALITY_NAMES[modality_name][i]
      if name not in dataset_names:
        if raise_exception:
          raise Exception(f"The dataset {name} is missing.")
        else:
          return False

  return True


def get_data(sonata_instance, dataset_name):
  dataset_path = f"{GROUP_0_PATH}/{dataset_name}"
  return sonata_instance[dataset_path][()]


def get_modality_data(sonata_instance, modality_name, downsampling = 1):
  if modality_name not in EXPECTED_MODALITY_NAMES:
    raise Exception(f"The modality {modality_name} does not exist.")

  data = {}
  for dataset_name in EXPECTED_MODALITY_NAMES[modality_name]:
    data[dataset_name] = get_data(sonata_instance, dataset_name)[::downsampling]
  
  return data


def modality_to_rab(sonata_instance, modality_name, output_file, downsampling = 1):
  # getting the data
  data = get_modality_data(sonata_instance, modality_name, downsampling)
  # building a Pandas Dataframe out of the data
  df = pd.DataFrame(data)
  # instanciating a RandomAccessBuffer object
  rab_instance = rab.RandomAccessBuffer()
  # adding an index to lookup
  rab_instance.addDataset("index", data = {
    "type": "CellRecordSeries",
    "dataset": modality_name,
  })
  # adding the dataframe to the RAB
  rab_instance.addDataset(modality_name, data = df, compress="gzip")
  rab_instance.write(output_file)


def euler_angle_to_rab(sonata_instance, output_file, downsampling = 1):
  """
  For the angles, they are encoded in the SONATA as quaternions, though in order to
  save a bit of space in the final RAB, we are converting them into euler angles
  """
  modality_name = "orientation"

  # getting the data
  data = get_modality_data(sonata_instance, modality_name, downsampling)
  
  # The quaternions
  q_w = data["orientation_w"]
  q_x = data["orientation_x"]
  q_y = data["orientation_y"]
  q_z = data["orientation_z"]

  # The eurler angles
  e_x = np.zeros_like(q_w)
  e_y = np.zeros_like(q_w)
  e_z = np.zeros_like(q_w)

  size = e_x.shape[0]
  
  for i in range(0, size):
    if np.isnan(q_w[i]):
      e_x[i] = np.nan
      e_y[i] = np.nan
      e_z[i] = np.nan
      continue

    r = Rotation.from_quat([q_x[i], q_y[i], q_z[i], q_w[i]])
    euler_angle = r.as_euler('xyz', degrees = False)
    e_x[i] = euler_angle[0]
    e_y[i] = euler_angle[1]
    e_z[i] = euler_angle[2]

  euler_data = {
    "orientation_x": e_x,
    "orientation_y": e_y,
    "orientation_z": e_z,
  }

  # building a Pandas Dataframe out of the data
  df = pd.DataFrame(euler_data)
  # instanciating a RandomAccessBuffer object
  rab_instance = rab.RandomAccessBuffer()
  # adding an index to lookup
  rab_instance.addDataset("index", data = {
    "type": "CellRecordSeries",
    "dataset": modality_name,
  })
  # adding the dataframe to the RAB
  rab_instance.addDataset(modality_name, data = df, compress="gzip")
  rab_instance.write(output_file)


def cell_type_to_rab(sonata_instance, output_file, downsampling = 1):
  """
  For the cell type, we need to export the metadata with the cell types
  """

  modality_name = "cell_type"
  data = get_modality_data(sonata_instance, modality_name, downsampling)
  # building a Pandas Dataframe out of the data
  df = pd.DataFrame(data)
  # instanciating a RandomAccessBuffer object
  rab_instance = rab.RandomAccessBuffer()
  # adding an index to lookup
  rab_instance.addDataset("index", data = {
    "type": "CellRecordSeries",
    "dataset": modality_name,
  },
  metadata = CELL_TYPE_ID_TO_CELL_TYPE_NAME)

  # adding the dataframe to the RAB
  rab_instance.addDataset(modality_name, data = df, compress="gzip")
  rab_instance.write(output_file)


def compute_stats(sonata_instance, hierarchy_filepath, output_file):
  cell_types = get_modality_data(sonata_instance, "cell_type")["cell_type"]
  regions = get_modality_data(sonata_instance, "region")["region_id"]

  unique_regions_in_sonata, count_cells_by_region = np.unique(regions, return_counts=True)
  unique_cell_types_in_sonata, count_cells_by_cell_type = np.unique(cell_types, return_counts=True)

  total_cell_count_by_cell_type = dict(zip(map(lambda i: CELL_TYPE_ID_TO_CELL_TYPE_NAME["lut"][i], unique_cell_types_in_sonata), count_cells_by_cell_type.tolist()))
  total_cell_count_by_leaf_region = dict(zip(unique_regions_in_sonata.tolist() , count_cells_by_region.tolist()))
  # above, leaf stands for "leaf as in the annotation volume"

  # print("cell_types", cell_types.shape)
  # print("regions", regions.shape)

  # where_cell_type_0 = np.where(cell_types == 0 )
  # print(len(where_cell_type_0[0]))
  cell_count_per_leaf_region_and_cell_type = {}

  nb_of_regions = len(unique_regions_in_sonata)
  counter = 0

  flag_value = len(CELL_TYPE_ID_TO_CELL_TYPE_NAME["lut"])

  for region_id in unique_regions_in_sonata:
    counter += 1
    # print(f"{counter}/{nb_of_regions} (region id: {region_id})")
    
    # this is an array that contains cell type ids only for the locations that are matchin the current region id.
    # The array is filled by default with a value that goes beyond the number of cell types
    # as a NO_DATA value so that we don't catch those in the next step
    mask_for_region = np.full_like(regions, fill_value = flag_value)
    where = regions == region_id
    mask_for_region[where] = cell_types[where]
    cell_types_only = np.delete(mask_for_region, mask_for_region == flag_value)
    unique, counts = np.unique(cell_types_only, return_counts=True)
    cell_count_per_leaf_region_and_cell_type[int(region_id)] = dict(zip(  map(lambda i: CELL_TYPE_ID_TO_CELL_TYPE_NAME["lut"][i], unique) , counts.tolist()))
  
  jsoncontent = json.loads(open(hierarchy_filepath, "r").read())

  # sometimes, the 1.json has its content in a "msg" sub prop (the original vesion has).
  # and some other versions don't. Here we deal with both
  if "msg" in jsoncontent:
    flat_tree = TreeIndexer.flattenTree(jsoncontent['msg'][0])
  else:
    flat_tree = TreeIndexer.flattenTree(jsoncontent)


  # For each region, we create a mask that contains all the sub regions
  region_counter = 0
  total_region = len(flat_tree)
  cell_count_per_region_and_cell_type = {} # for non leaf regions
  for region_id in flat_tree:

    region_counter += 1
    region_node = flat_tree[region_id]

    # print(f"{region_counter}/{total_region} - [{region_id}] {flat_tree[region_id]['name']}")
    cell_count_by_cell_type_for_this_region = {}
    
    descendant_including_current = region_node["_descendants"][:]
    descendant_including_current.append(region_id)

    total_for_this_region = 0
    for sub_region_id in descendant_including_current:
      if sub_region_id not in cell_count_per_leaf_region_and_cell_type:
        continue

      for cell_type_name in cell_count_per_leaf_region_and_cell_type[sub_region_id]:
        if cell_type_name not in cell_count_by_cell_type_for_this_region:
          cell_count_by_cell_type_for_this_region[cell_type_name] = 0
        cell_count_by_cell_type_for_this_region[cell_type_name] += cell_count_per_leaf_region_and_cell_type[sub_region_id][cell_type_name]
        total_for_this_region += cell_count_per_leaf_region_and_cell_type[sub_region_id][cell_type_name]

    # if total_for_this_region:
    cell_count_per_region_and_cell_type[region_id] = {
      "id": region_id,
      "total": total_for_this_region,
      "ratioToTotal": float(total_for_this_region / len(regions)),
      "perCellType": cell_count_by_cell_type_for_this_region,
    }

    ratio_to_total_cell_type = {}
    for cell_type in cell_count_by_cell_type_for_this_region:
      ratio_to_total_cell_type[cell_type] = float(cell_count_by_cell_type_for_this_region[cell_type] / total_cell_count_by_cell_type[cell_type])

    cell_count_per_region_and_cell_type[region_id]["ratioToTotaPerCellType"] = ratio_to_total_cell_type

  f = open(output_file, "w+")
  f.write(json.dumps(cell_count_per_region_and_cell_type, indent = 2))
  f.close()



def main():
  args = parse_args(sys.argv[1:])

  sonata = h5py.File(args.sonata_filepath,"r")
  
  if not check_sonata(sonata, raise_exception = False):
    print("This SONATA is invalid (properties are missing")
    exit(1)

  print("Computing statistics...")
  compute_stats(sonata, hierarchy_filepath = args.hierarchy_filepath, output_file = args.out_metadata)

  print("Exporting positions...")
  modality_to_rab(sonata, "position", args.out_position, downsampling = args.downsampling_factor)

  print("Exporting orientations...")
  euler_angle_to_rab(sonata, args.out_orientation, downsampling = args.downsampling_factor)

  print("Exporting regions...")
  modality_to_rab(sonata, "region", args.out_region, downsampling = args.downsampling_factor)

  print("Exporting cell types...")
  cell_type_to_rab(sonata, args.out_cell_type, downsampling = args.downsampling_factor)


if __name__ == "__main__":
  main()

