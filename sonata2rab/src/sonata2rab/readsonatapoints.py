import h5py
import numpy as np
import randomaccessbuffer as rab
import pandas as pd



GROUP_0_PATH = "nodes/atlas_cells/0"
EXPECTED_MODALITY_NAMES = {
    "position": ["x", "y", "z"],
    "orientation": ["orientation_w", "orientation_x", "orientation_y", "orientation_z"],
    "region": ["region_id"],
    "cell_type": ["cell_type"],
  }


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
    data[dataset_name] = np.nan_to_num(get_data(sonata_instance, dataset_name)[::downsampling])
  
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
    "modality": modality_name,
  })
  # adding the dataframe to the RAB
  rab_instance.addDataset(modality_name, data = df, compress="gzip")
  rab_instance.write(output_file)


if __name__ == "__main__":
    sonata = h5py.File("/Users/lurie/Documents/BBP/projects/blue_brain_atlas_website_exporter/sonata2rab/sample_data/cell_positions_hybrid.h5","r")
    
    if not check_sonata(sonata, raise_exception = False):
      print("This SONATA ain't no good.")
      exit()

    modality_to_rab(sonata, "position", "/Users/lurie/Documents/BBP/projects/blue_brain_atlas_website_exporter/sonata2rab/sample_data/position.rab")
    modality_to_rab(sonata, "orientation", "/Users/lurie/Documents/BBP/projects/blue_brain_atlas_website_exporter/sonata2rab/sample_data/orientation.rab")
    modality_to_rab(sonata, "region", "/Users/lurie/Documents/BBP/projects/blue_brain_atlas_website_exporter/sonata2rab/sample_data/region.rab")
    modality_to_rab(sonata, "cell_type", "/Users/lurie/Documents/BBP/projects/blue_brain_atlas_website_exporter/sonata2rab/sample_data/cell_type.rab")
