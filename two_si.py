#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sirius, json
import numpy as np
from sirius_md import cpmd_verlet, verlet


# ## Setup the unit cell

# In[7]:


lat = 10.26* np.eye(3)
gk_cutoff = 6
pw_cutoff = 12
inp={
    "parameters" : {
        "xc_functionals" : ["XC_LDA_X", "XC_LDA_C_PZ"],
        "electronic_structure_method" : "pseudopotential",
        "pw_cutoff" : pw_cutoff,
        "gk_cutoff" : gk_cutoff,
        "num_fv_states": 4,
        "use_symmetry": True,
        "gamma_point": False
    },
    "control" : {
        "verbosity" : 0
    }
}
ctx = sirius.Simulation_context(json.dumps(inp))
ctx.unit_cell().set_lattice_vectors(*lat)
ctx.unit_cell().add_atom_type('Si','pseudos/Si.json')
ctx.unit_cell().add_atom('Si', [0,0,0])
ctx.unit_cell().add_atom('Si', [0.25,]*3)
ctx.initialize()


# In[8]:


kset = sirius.K_point_set(ctx, [1,1,1], [0,0,0], True)


# In[9]:


dft = sirius.DFT_ground_state(kset)
dft.initial_state()
tol = 1e-9
dft.find(tol, tol, tol, 1000, False)


# In[12]:


x0 = sirius.atom_positions(ctx.unit_cell())
v0 = np.zeros_like(x0)

input_vars = {'parameters': {'dt': 4, 'N': 1000, 'me': 300, 'T': 0}}

for step_data in cpmd_verlet.cpmd_verlet_raw(input_vars, kset):
    print(step_data)
print('end')

