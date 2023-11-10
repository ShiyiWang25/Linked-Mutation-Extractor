
import sys

from . import lime, utils

args = utils._parse_args(sys.argv[1:])

lime.lime(args)
