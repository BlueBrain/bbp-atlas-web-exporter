"""
Get the actual version of the package as it is installed
that depends on the git tag. Gives the same as 'pip show'.

(As opposed as a hardcoded version number that the 
developer has to maintain)
"""

import pkg_resources
# This is the name of your current package itself
from blue_brain_atlas_web_exporter import __name__
__version__ = pkg_resources.get_distribution(__name__).version
