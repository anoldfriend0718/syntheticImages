
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

with open(parameters["rock_simulation"]["simulation_rock_file"], "r") as fp:
    raw_data = np.fromfile(fp, dtype=np.uint8)
simulation_rock = raw_data.reshape(
    parameters["rock_simulation"]["simulation_rock_nx"],
    parameters["rock_simulation"]["simulation_rock_ny"])
generator = SedimentBearingRockGenerator()
sediment_bearing_rock = generator.deposit(simulation_rock, parameters)
generator.write_in_openfoam_format(sediment_bearing_rock, parameters)
