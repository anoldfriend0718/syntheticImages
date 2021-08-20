
import json
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors

sys.path.append(f"/home/anoldfriend/Workspace/MyRepo/hbs_syntenic_images/syntheticImages/pyScripts")
from sedimentRock import SedimentBearingRockGenerator

with open("./parameters.json", "r") as fp:
    parameters = json.load(fp)
generator = SedimentBearingRockGenerator()
simulation_rock = generator.simulate(parameters)

sediment_bearing_rock = generator.deposit(simulation_rock, parameters)
generator.write_in_openfoam_format(sediment_bearing_rock, parameters)
