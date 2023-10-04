#!/usr/bin/env python
# coding: utf-8

import sirius, json
from sirius_md import cpmd_verlet
from sirius.coefficient_array import zeros_like
import numpy as np

# ## Setup the unit cell
# In[3]:
lat = 10.26 * np.eye(3)
gk_cutoff = 3
pw_cutoff = 6
inp = {
    "parameters": {
        "xc_functionals": ["XC_LDA_X", "XC_LDA_C_PZ"],
        "electronic_structure_method": "pseudopotential",
        "pw_cutoff": pw_cutoff,
        "gk_cutoff": gk_cutoff,
    },
    "control": {"verbosity": 2},
}
ctx = sirius.Simulation_context(json.dumps(inp))
ctx.unit_cell().set_lattice_vectors(*lat)
ctx.unit_cell().add_atom_type("Si", "pseudos/Si.json")
ctx.unit_cell().add_atom("Si", [0, 0, 0])
ctx.unit_cell().add_atom(
    "Si",
    [
        0.2532,
    ]
    * 3,
)
ctx.initialize()

kset = sirius.K_point_set(ctx, [1, 1, 1], [0, 0, 0], True)

dft = sirius.DFT_ground_state(kset)
tol = 1e-9
dft.find(tol, tol, tol, 1000, False)

x0 = sirius.atom_positions(ctx.unit_cell())
v0 = np.zeros_like(x0)

input_vars = {'parameters': {'dt': 10, 'N': 1000, 'me': 300, 'T': 0}}

for step_data in cpmd_verlet.cpmd_verlet_raw(input_vars, kset):
    print(step_data)
print('end')
