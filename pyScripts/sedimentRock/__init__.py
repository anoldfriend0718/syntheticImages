import sys
import os.path as path
sys.path.append(path.dirname(path.abspath(__file__))) #add sys path, avoiding  modules in this package cannot find their peer module

from .sediment_bearing_rock_generator import SedimentBearingRockGenerator
