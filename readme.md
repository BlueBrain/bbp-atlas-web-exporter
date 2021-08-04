# BlueBrainAtlas WebExporter

The web exporter is a suite of tools to facilitate the conversion and export of the Blue Brain Atlas pipeline to web-friendly datasets to be consumed by the Atlas web app.  

Once installed, this package exposes two CLIs:
- `cellexport`
- `parcellationexport`

# Installation
`pip install git+https://bbpgitlab.epfl.ch/dke/apps/blue_brain_atlas_web_exporter.git`  
This will automatically install the CLIs listed below and make them available in your `PATH`.

# Usage
## cellexport
The CLI `cellexport` converts a cell-position SONATA file into a set of RandomAccessBuffer files (RAB) as well as a JSON metadata file containing a few computed metrics.  
Here is the list of arguments:
- `--help` displays the help menu
- `--version` displays the version number
- `--hierarchy <FILE PATH>` input path to the brain region hierarchy file (commonly called *1.json*) — **mandatory**
- `--sonata <FILE PATH>` input path to the SONATA file containing cell positions/orientation/brain-regions/celltypes — **mandatory**
- `--downsampling-factor <NUMBER>` factor to downsample the data — **optional** (default: `1`, no downsampling. Example: a factor of 100 will keep only one value every 100)
- `--out-postion <FILE PATH>` output path to the RAB file that contains the cell positions — **mandatory**
- `--out-cell-type <FILE PATH>` output path to the RAB file that contains the cell types — **mandatory**
- `--out-orientation <FILE PATH>` output path to the RAB file that contains the cell orientations — **mandatory**
- `--out-region <FILE PATH>` output path to the RAB file that contains the cell region IDs — **mandatory**
- `--out-metadata <FILE PATH>` output path the JSON file that contains the metadata — **mandatory**


## parcellationexport
The CLI `parcellationexport` computes all the meshes (`.OBJ` files) as well as the binary masks (`.nrrd` files) for all the regions available in both the hierarchy file (input) and the annotation volume (input). Many intermediate regions from the hierarchy file do not have any voxel representation and instead are composed of voxels from their child region, thus, this CLI computes all the aggregations necessary to recreate binary masks for those regions. This CLI also output a JSON metadata file with metrics regarding the volume and volume ration for all the regions.  
Here is the list of arguments:  
- `--help` displays the help menu
- `--version` displays the version number
- `--hierarchy <FILE PATH>` input path to the brain region hierarchy file (commonly called *1.json*) — **mandatory**
- `--parcellation-volume <FILE PATH>` input path to the parcellation file format (*.nrrd*) — **mandatory**
- `--out-mesh-dir <DIR PATH>` output path for the directory where the region/parcellation mesh files (*.obj*) will be written. Each file will have the name corresponding to the ID of the region (example: `997.obj` for the whole brain) — **mandatory**
- `--out-mask-dir <DIR PATH>` output path for the directory where the region/parcellation binary mask volumes (*.nrrd*). Each file will have the name corresponding to the ID of the region (example: `997.nrrd` for the whole brain) — **mandatory**
- `--out-metadata <FILE PATH>` output path to the JSON file that contains the metadata - **mandatory**

# Contributors
- Jonathan Lurie <jonathan.lurie@epfl.ch> (DKE team)


# License
TODO