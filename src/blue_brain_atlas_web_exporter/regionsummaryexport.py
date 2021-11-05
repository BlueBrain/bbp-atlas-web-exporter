# -*- coding: utf-8 -*-
"""

"""

import argparse
import sys
import json
import os
import platform
import subprocess
import re
import nrrd
import numpy as np
import blue_brain_atlas_web_exporter.TreeIndexer as TreeIndexer
import blue_brain_atlas_web_exporter
from blue_brain_atlas_web_exporter import __version__

def writeMetadata(data, filepath):
    metadata_file = open(filepath, 'w')
    metadata_file.write(json.dumps(data, ensure_ascii = False, indent = 2))
    metadata_file.close()

def parse_args(args):
    """Parse command line parameters

    Args:
      args ([str]): command line parameters as list of strings

    Returns:
      :obj:`argparse.Namespace`: command line parameters namespace
    """
    parser = argparse.ArgumentParser(
        description="Aggregates region informations and export a summary.")
    parser.add_argument(
        "--version",
        action="version",
        version="regionsummaryexport {ver}".format(ver=__version__))

    parser.add_argument(
        "--hierarchy",
        dest="hierarchy",
        required=True,
        metavar="<FILE PATH>",
        help="The hierarchy JSON file, sometimes called 1.json")

    parser.add_argument(
        "--cell-metadata",
        dest="cell_metadata",
        required=True,
        metavar="<FILE PATH>",
        help="The JSON file containing the cell and cell density metadata for all the regions")

    parser.add_argument(
        "--parcellation-metadata",
        dest="parcellation_metadata",
        required=True,
        metavar="<FILE PATH>",
        help="The JSON file containing the parcellation informations such as volume and neighbours")

    parser.add_argument(
        "--link-metadata",
        dest="link_metadata",
        required=True,
        metavar="<FILE PATH>",
        help="The JSON file containing the links between regions and their resources in Nexus (@id of meshes, masks, atlasRelease, etc.")

    parser.add_argument(
        "--region-prefix",
        dest="region_prefix",
        default="mba:",
        help="The prefix added to the region IDs (default 'mba:')")

    parser.add_argument(
        "--resource-type",
        dest="resource_type",
        nargs='+',
        default="RegionSummary",
        help="The @type for such resources. Whitespaces separted if many (default: 'RegionSummary')")

    parser.add_argument(
        "--output",
        dest="out_metadata",
        required=True,
        metavar="<FILE PATH>",
        help="Path to the output metadata file (json)")


    return parser.parse_args(args)


def main():
    """Main entry point allowing external calls

    Args:
      args ([str]): command line parameter list
    """
    args = parse_args(sys.argv[1:])

    hierarchy_filepath = args.hierarchy
    cell_metadata_filepath = args.cell_metadata
    parcellation_metadata_filepath = args.parcellation_metadata
    link_metadata_filepath = args.link_metadata
    out_metadata_filepath = args.out_metadata
    region_prefix = args.region_prefix
    resource_type = args.resource_type

    # loading hierarchy json
    hierarchy = json.loads(open(hierarchy_filepath, "r").read())

    # sometimes, the 1.json has its content in a "msg" sub prop (the original vesion has).
    # and some other versions don't. Here we deal with both
    if "msg" in hierarchy:
        flat_tree = TreeIndexer.flattenTree(hierarchy['msg'][0])
    else:
        flat_tree = TreeIndexer.flattenTree(hierarchy)

    # Reading the other metadata input files
    cell_metadata = json.loads(open(cell_metadata_filepath, "r").read())
    parcellation_metadata = json.loads(open(parcellation_metadata_filepath, "r").read())
    link_metadata = json.loads(open(link_metadata_filepath, "r").read())

    all_keys = set()
    all_keys.update(cell_metadata.keys())
    all_keys.update(parcellation_metadata.keys())
    all_keys.update(parcellation_metadata.keys())

    summary = {}

    for region_id in all_keys:
      present = True

      

      if int(region_id) not in flat_tree:
        present = False
        # print(f"⚠️ region {region_id} is missing from the hierarchy file")

      if region_id not in cell_metadata:
        present = False
        # print(f"⚠️ region {region_id} is missing from the cell metadata file")

      if region_id not in parcellation_metadata:
        present = False
        # print(f"⚠️ region {region_id} is missing from the parcellation metadata file")

      if region_id not in link_metadata:
        present = False
        # print(f"⚠️ region {region_id} is missing from the link metadata file")

      if not present:
        continue

      hierarchy_node = flat_tree[int(region_id)]

      adj_to = parcellation_metadata[region_id]["adjacentTo"]
      cont_w = parcellation_metadata[region_id]["continuousWith"]

      cell_meta = cell_metadata[region_id]
      cell_counts = cell_meta["perCellType"]
      # add an entry for the nsg:Cell total count
      cell_counts["nsg:Cell"] = cell_meta["total"]
      
      cell_ratios = cell_meta["ratioToTotaPerCellType"]
      # add an entry for the nsg:Cell total ratio
      cell_ratios["nsg:Cell"] = cell_meta["ratioToTotal"]

      cur_reg_dict = {
        # data from ontology
        "@type": resource_type,
        "name": hierarchy_node["name"],
        "aconym": hierarchy_node["acronym"],
        "color": hierarchy_node["color_hex_triplet"],
        "description": f"Summary of informations regarding the region {hierarchy_node['name']} as seen in the atlas release {link_metadata[region_id]['atlasRelease']['@id']}.",
        
        # AtlasRelease
        "atlasRelease": {
          "@id": link_metadata[region_id]['atlasRelease']['@id']
        },

        # SRS and prefixed region id
        "brainLocation": {
          "atlasSpatialReferenceSystem": {
            "@id": link_metadata[region_id]['atlasSpatialReferenceSystem']['@id'],
            "@type": [
              "AtlasSpatialReferenceSystem",
              "BrainAtlasSpatialReferenceSystem"
            ]
          },
          "brainRegion": {
            "@id": f"{region_prefix}{region_id}"
          }
        },

        # Link to the mesh resource
        "mesh": {
          "@id": link_metadata[region_id]['mesh']['@id'],
        },

        # link to mask resource
        "mask": {
          "@id": link_metadata[region_id]['mask']['@id'],
        },

        # volume information (in SI units)
        "volume": {
          "total": {
            "size": parcellation_metadata[region_id]["regionVolume"],
            "unitCode": "cubic micrometer"
          },

          "ratio": parcellation_metadata[region_id]["regionVolumeRatioToWholeBrain"]
        },


        # regions that are sharing a border with this current one
        # (regardless of a possible common function)
        "adjacentTo": list(map(lambda k:  { "@id": f"{region_prefix}{k}", "ratio": adj_to[k]} , adj_to)),


        # Subset of adjacentTo but only populated with regions that, in addition to sharing a border,
        # share a comon function, is a continuation within the same region.
        "continuousWith": [
          list(map(lambda el:  { "@id": f"{region_prefix}{el}" } , cont_w)),
        ],

        "layers": parcellation_metadata[region_id]["layers"],

        # cell_meta = cell_metadata[region_id]
        # cell_counts = cell_meta["perCellType"]
        # cell_ratios = cell_meta["ratioToTotaPerCellType"]

        "cellMetrics": {
          "cellCount": list(map(lambda k:  { "cellType": k, "value": cell_counts[k]} , cell_counts)),

          # "cellCount": [
          #   {
          #     "cellType": "nsg:Cell", # This is actually about a metric for all the cell,
          #     "value": 2918225,       # the specialized cell types come below.
          #   },
          #   {
          #     "cellType": "nsg:ExcitatoryNeuron",
          #     "value": 304919
          #   },
          #   {
          #     "cellType": "nsg:InhibitoryNeuron",
          #     "value": 1487618,
          #   },
          #   {
          #     "cellType": "nsg:Microglia",
          #     "value": 7835
          #   },
          #   {
          #     "cellType": "nsg:Astrocyte",
          #     "value": 398519
          #   },
          #   {
          #     "cellType": "nsg:Oligodendrocyte",
          #     "value": 719334
          #   }
          # ],

          # ratios are with the whole brain cell count as reference (in this particular atlas)
          "cellRatio": list(map(lambda k:  { "cellType": k, "value": cell_ratios[k]} , cell_ratios)),
          # [
          #   {
          #     "cellType": "nsg:Cell",
          #     "@value": 0.02621629702531575
          #   },
          #   {
          #     "cellType": "nsg:ExcitatoryNeuron",
          #     "@value": 0.03833953509074801
          #   },
          #   {
          #     "cellType": "nsg:InhibitoryNeuron",
          #     "@value": 0.022519446251155515
          #   },
          #   {
          #     "cellType": "nsg:Microglia",
          #     "@value": 0.001374382733753384
          #   },
          #   {
          #     "cellType": "nsg:Astrocyte",
          #     "@value": 0.027730323803069116
          #   },
          #   {
          #     "cellType": "nsg:Oligodendrocyte",
          #     "@value": 0.04175126707626208
          #   }
          # ]
        }


      }

      summary[region_id] = cur_reg_dict

    writeMetadata(summary, out_metadata_filepath)

    

if __name__ == "__main__":
    main()
