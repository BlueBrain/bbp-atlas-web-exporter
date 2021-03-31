import h5py
import numpy as np


GROUP_0_PATH = "nodes/atlas_cells/0"
EXPECTED_DATASET_NAMES = {
    "position": ["x", "y", "z"],
    "orientation": ["orientation_w", "orientation_x", "orientation_y", "orientation_z"],
    "region": ["region_id"],
    "cell_type": ["cell_type"],
  }


# Check that all the expected datasets are present in the SONATA file
def check_sonata(sonata_instance, raise_exception = True):
  group_0 = sonata_instance[GROUP_0_PATH]
  dataset_names = list(group_0.keys())

  for category_name in EXPECTED_DATASET_NAMES:
    for i in range(0, len(EXPECTED_DATASET_NAMES[category_name])):
      name = EXPECTED_DATASET_NAMES[category_name][i]
      if name not in dataset_names:
        if raise_exception:
          raise Exception(f"The dataset {name} is missing.")
        else:
          return False

  return True


def get_data(sonata_instance, dataset_name):
  dataset_path = f"{GROUP_0_PATH}/{dataset_name}"
  return sonata_instance[dataset_path].value


def get_position(sonata_instance, downsampling = 1):
  return {
    "x": get_data(sonata_instance, "x")[::downsampling],
    "y": get_data(sonata_instance, "y")[::downsampling],
    "z": get_data(sonata_instance, "z")[::downsampling],
  }
  

# also replaces NaNs with zeros
def get_orientation(sonata_instance, downsampling = 1):
  return {
    "orientation_x": np.nan_to_num(get_data(sonata_instance, "orientation_x")[::downsampling]),
    "orientation_y": np.nan_to_num(get_data(sonata_instance, "orientation_y")[::downsampling]),
    "orientation_z": np.nan_to_num(get_data(sonata_instance, "orientation_z")[::downsampling]),
  }


def get_region_id(sonata_instance, downsampling = 1):
  return {
    "region_id": get_data(sonata_instance, "region_id")[::downsampling],
  }


def get_cell_type(sonata_instance, downsampling = 1):
  return {
    "cell_type": get_data(sonata_instance, "cell_type")[::downsampling],
  }


if __name__ == "__main__":
    sonata = h5py.File("/Users/lurie/Documents/BBP/projects/blue_brain_atlas_website_exporter/sonata2rab/sample_data/cell_positions_hybrid.h5","r")
    
    if not check_sonata(sonata, raise_exception = False):
      print("This SONATA ain't no good.")
      exit()

    position = get_position(sonata)
    orientation = get_orientation(sonata)
    region_id = get_region_id(sonata)
    cell_type = get_cell_type(sonata)
    
    print(position)
    print(orientation)
    print(region_id)
    print(cell_type)


# print( group_0["x"] )

# xpos = np.float32(h5file["x"].value)