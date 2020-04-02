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
import os


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
        "--in-dir",
        dest="in_dir",
        required=True,
        metavar="<DIRECTORY PATH>",
        help="The input directory with OBJs to smooth")

    parser.add_argument(
        "--out-dir",
        dest="out_dir",
        required=True,
        metavar="<DIRECTORY PATH>",
        help="The output directory to create smooth OBJs")


    return parser.parse_args(args)



def main():
    """Main entry point allowing external calls

    Args:
      args ([str]): command line parameter list
    """
    args = parse_args(sys.argv[1:])
    out_dir = args.out_dir
    in_dir = args.in_dir

    # create out_dir if inexistant
    try:
        os.makedirs(out_dir)
    except FileExistsError as e:
        pass


    blender_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "blender_bin", "blender-2.79-macOS-10.6", "blender.app", "Contents", "MacOS", "blender")
    script_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "smoother_plugin_script.py")

    command = blender_path + " -b -P " + script_path + " " + in_dir + " " + out_dir
    print(command)
    os.system(command)

if __name__ == "__main__":
    main()
