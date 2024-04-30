import json
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

case = 'replicate_1d'
n = 200
input_path = '/Users/ashtoncole/Documents/School/COE_379L_parallel/Parallel_Computing/final/test'
output_path = '/Users/ashtoncole/Documents/School/COE_379L_parallel/Parallel_Computing/final/doc/images'

# Metadata
file_meta = open(f'{input_path}/{case}/metadata.json')
metadata = json.load(file_meta)
file_meta.close()

# Case data
file_case = open(f'{input_path}/{case}/{n}/t')
t = float(file_case.readline())
file_case.close()

file_case = open(f'{input_path}/{case}/{n}/z')
file_case.readline()
file_case.readline()
z = file_case.read().split('\n')
z = np.reshape(
    np.array([float(val) for val in z[0:(len(z) - 1)]]),
    (metadata["SpaceDomain"]["nx"], metadata["SpaceDomain"]["ny"], 3)
)
file_case.close()

rho_p = z[:, :, 0]
u_p = z[:, :, 1]
v_p = z[:, :, 2]

# Plot
xs = np.linspace(
    metadata["SpaceDomain"]["xa"],
    metadata["SpaceDomain"]["xb"],
    metadata["SpaceDomain"]["nx"]
)
ys = np.linspace(
    metadata["SpaceDomain"]["ya"],
    metadata["SpaceDomain"]["yb"],
    metadata["SpaceDomain"]["ny"]
)
print(ys)
Y, X = np.meshgrid(ys, xs)

fig, ax = plt.subplots()

# rho_p
ax.pcolor(X, Y, rho_p)
fig.savefig(f'{output_path}/{case}_rho_p_{n}.png')

# u_p
ax.pcolor(X, Y, u_p)
fig.savefig(f'{output_path}/{case}_u_p_{n}.png')

# v_p
ax.pcolor(X, Y, v_p)
fig.savefig(f'{output_path}/{case}_v_p_{n}.png')

# Save