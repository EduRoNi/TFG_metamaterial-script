'''
Based on the data extracted from the other scripts, the yield surface f(sigma) is established as a 
failure and/or plastic flow criterion for the new material defined by homogenizing the metamaterial structure.
'''
import numpy as np
import matplotlib.pyplot as plt

import plotly.graph_objects as go
from skimage import measure
import sys

import shutil
import os
import pandas as pd

######################################################################
##                             SETTING                              ##
######################################################################
# Fonts and style
plt.rcParams.update({
    "text.usetex": False,  # Cambiar a True si tienes LaTeX instalado en el sistema
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman", "Times New Roman", "DejaVu Serif"],
    "font.size": 12,
    "axes.edgecolor": "black",
    "axes.linewidth": 0.8,
    "xtick.direction": "in",
    "ytick.direction": "in",
}) 
######################################################################
##                              DATA                                ##
######################################################################
# Read file
if len(sys.argv) >= 3:
    file_name = sys.argv[1]
    reports_dir = sys.argv[2]
    temporal = sys.argv[3]

reports_dir_cell = os.path.join(reports_dir, file_name)
file_path = os.path.join(reports_dir_cell, f"reporte_{file_name}.xlsx")

meta = pd.read_excel(file_path, sheet_name="Info")
Sigma_x = meta.loc[meta["Parámetro"] == "Límite elast_X (MPa)", "Valor"].values[0] # [MPa]
Sigma_y = meta.loc[meta["Parámetro"] == "Límite elast_Y (MPa)", "Valor"].values[0] # [MPa]
Sigma_z = meta.loc[meta["Parámetro"] == "Límite elast_Z (MPa)", "Valor"].values[0] # [MPa]
Sigma_xy = meta.loc[meta["Parámetro"] == "Límite elast_YX (MPa)", "Valor"].values[0] # [MPa]
Sigma_xz = meta.loc[meta["Parámetro"] == "Límite elast_XZ (MPa)", "Valor"].values[0] # [MPa]
Sigma_yz = meta.loc[meta["Parámetro"] == "Límite elast_YZ (MPa)", "Valor"].values[0] # [MPa]

datos = np.array([Sigma_x, Sigma_y, Sigma_z, Sigma_xy, Sigma_xz, Sigma_yz])

# YIELD STRESS
# Continuum body Ti-6Al-4V
SY = 945 # MPa

######################################################################
##                             SURFACE                              ##
######################################################################
mat_real = 0 # 0-Do not plot the yield surface of the real material. 1-Plot the yield surface of the real material.
# f(sigma) = F*(s22 - s33)**2 + G*(s33 - s11)**2 + H*(s11 - s22)**2 + 2*L*s23**2 + 2*M*s13**2 + 2*N*s12**2 = 1
# CONSTANTES
F = 0.5*( 1/Sigma_y**2 + 1/Sigma_z**2 - 1/Sigma_x**2)
G = 0.5*( 1/Sigma_z**2 + 1/Sigma_x**2 - 1/Sigma_y**2)
H = 0.5*( 1/Sigma_x**2 + 1/Sigma_y**2 - 1/Sigma_z**2)
R = 0.5*(1/SY**2) # for real continuum body

N = 1/(2*Sigma_xy**2)
M = 1/(2*Sigma_xz**2)
L = 1/(2*Sigma_yz**2)

print("[", end='')
print(F, G, H, N, M, L, sep=',', end='')
print("]")

limit = max(datos)*1.25
if mat_real == 1:
    limit = SY*1.25 # max(Sigma_x, Sigma_y, Sigma_z)*1.25
s_puntos = np.linspace(-limit, limit, 200)

# ISOLINES 2D
fig, ax1 = plt.subplots(figsize=(7, 6))
# s3=0
s1, s2= np.meshgrid(s_puntos, s_puntos)
f = F*(s2)**2 + G*(s1)**2 + H*(s1 - s2)**2 - 1
cs = ax1.contour(s1, s2, f, levels=[0], colors='blue', linewidth=2)
ax1.set_title('Principle Plane 1-2', fontsize=18)
ax1.set_xlabel('$\sigma_1$ [MPa]', fontsize=14, labelpad=8)
ax1.set_ylabel('$\sigma_2$ [MPa]', fontsize=14, labelpad=8)
ax1.axis('equal')
ax1.grid(True, linestyle='-', linewidth=0.5, color='lightgray')
artists, _ = cs.legend_elements()

if mat_real == 1:
    g = R*(s2)**2 + R*(s1)**2 + R*(s1 - s2)**2 - 1 
    gs = ax1.contour(s1, s2, g, levels=[0], colors='red')
    artistss, _ = gs.legend_elements()
    ax1.legend([artists[0], artistss[0]], ['Metamaterial', 'Bulk Solid'], 
           loc='upper left', bbox_to_anchor=(0.04, 0.96), edgecolor='black', framealpha=1)
elif mat_real ==0 :
    ax1.legend([artists[0]], ['Metamaterial'], 
           loc='upper left', bbox_to_anchor=(0.04, 0.96), edgecolor='black', framealpha=1)

plt.tight_layout()
plt.savefig(f"Plano_principal_1-2.png", dpi=300)
plt.close()
shutil.move(f"Plano_principal_1-2.png", os.path.join(reports_dir_cell, f"Plano_principal_1-2.png")) # move to correct directory


fig, ax2 = plt.subplots(figsize=(7, 6))
# s2=0
s1, s3= np.meshgrid(s_puntos, s_puntos)
f = F*(s3)**2 + G*(s1-s3)**2 + H*(s1)**2 - 1 

cs = ax2.contour(s1, s3, f, levels=[0], colors='blue', linewidth=2)
ax2.set_title('Principle Plane 1-3', fontsize=18)
ax2.set_xlabel('$\sigma_1$ [MPa]', fontsize=14, labelpad=8)
ax2.set_ylabel('$\sigma_3$ [MPa]', fontsize=14, labelpad=8)
ax2.axis('equal')
ax2.grid(True, linestyle='-', linewidth=0.5, color='lightgray')
artists, _ = cs.legend_elements()

if mat_real == 1:
    g = R*(s3)**2 + R*(s1-s3)**2 + R*(s1)**2 - 1
    gs = ax2.contour(s1, s3, g, levels=[0], colors='red')
    artistss, _ = gs.legend_elements()
    ax2.legend([artists[0], artistss[0]], ['Metamaterial', 'Bulk Solid'], 
           loc='upper left', bbox_to_anchor=(0.04, 0.96), edgecolor='black', framealpha=1)
elif mat_real == 0:
    ax2.legend([artists[0]], ['Metamaterial'], 
           loc='upper left', bbox_to_anchor=(0.04, 0.96), edgecolor='black', framealpha=1)






plt.tight_layout()
plt.savefig(f"Plano_principal_1-3.png", dpi=300)
plt.close()
shutil.move(f"Plano_principal_1-3.png", os.path.join(reports_dir_cell, f"Plano_principal_1-3.png")) # move to correct directory


fig, ax3 = plt.subplots(figsize=(7, 6))
# s1=0
s2, s3= np.meshgrid(s_puntos, s_puntos)

cs = ax3.contour(s2, s3, f, levels=[0], colors='blue', linewidth=2)
ax3.set_title('Principle Plane 2-3', fontsize=18)
ax3.set_xlabel('$\sigma_2$ [MPa]', fontsize=14, labelpad=8)
ax3.set_ylabel('$\sigma_3$ [MPa]', fontsize=14, labelpad=8)
ax3.axis('equal')
ax3.grid(True, linestyle='-', linewidth=0.5, color='lightgray')
artists, _ = cs.legend_elements()

f = F*(s2-s3)**2 + G*(s3)**2 + H*(s2)**2 - 1 
if mat_real == 1:
    g = R*(s2-s3)**2 + R*(s3)**2 + R*(s2)**2 - 1 
    gs = ax3.contour(s2, s3, g, levels=[0], colors='red')
    artistss, _ = gs.legend_elements()
    ax3.legend([artists[0], artistss[0]], ['Metamaterial', 'Bulk Solid'], 
           loc='upper left', bbox_to_anchor=(0.04, 0.96), edgecolor='black', framealpha=1)
elif mat_real == 0:
    ax3.legend([artists[0]], ['Metamaterial'], 
           loc='upper left', bbox_to_anchor=(0.04, 0.96), edgecolor='black', framealpha=1)

plt.tight_layout()
plt.savefig(f"Plano_principal_2-3.png", dpi=300)
plt.close()
shutil.move(f"Plano_principal_2-3.png", os.path.join(reports_dir_cell, f"Plano_principal_2-3.png")) # move to correct directory


fig, ax4 = plt.subplots(figsize=(7, 6))
# plane pi
U, V = np.meshgrid(s_puntos, s_puntos)
# Inverse transformation: from Pi plane (u,v) to stress space (s1, s2, s3)
# Assuming zero hydrostatic component for the projection
s1_pi = -U/np.sqrt(2) - V/np.sqrt(6)
s2_pi =  U/np.sqrt(2) - V/np.sqrt(6)
s3_pi =  2*V/np.sqrt(6)
# Hill equation
f_pi = F*(s2_pi - s3_pi)**2 + G*(s3_pi - s1_pi)**2 + H*(s1_pi - s2_pi)**2 - 1

cs = ax4.contour(U, V, f_pi, levels=[0], colors='blue')
ax4.set_title('Plane $\pi$', fontsize=18)
ax4.set_xlabel('$u$', fontsize=14, labelpad=8)
ax4.set_ylabel('$u$', fontsize=14, labelpad=8)
ax4.axis('equal')
ax4.grid(True, linestyle='-', linewidth=0.5, color='lightgray')

if mat_real == 1:
    g_pi = R*(s2_pi - s3_pi)**2 + R*(s3_pi - s1_pi)**2 + R*(s1_pi - s2_pi)**2 - 1
    gs = ax4.contour(U, V, g_pi, levels=[0], colors='red')


# Axis s1, s2, s3 for reference
for angle in [210, 330, 90]: # Angles de s1, s2, s3 proyected
    rad = np.deg2rad(angle)
    ax4.plot([0, limit*np.cos(rad)], [0, limit*np.sin(rad)], color="black", alpha=0.5)

plt.tight_layout()
plt.savefig(f"Plano_pi.png", dpi=300)
plt.close()
shutil.move(f"Plano_pi.png", os.path.join(reports_dir_cell, f"Plano_pi.png")) # move to correct directory

# SURFACE 3D
# Grid
lim = limit
res = 20 
s1, s2, s3 = np.meshgrid(np.linspace(-lim, lim, res),
                         np.linspace(-lim, lim, res),
                         np.linspace(-lim, lim, res), indexing='ij')


# Hill equation for principle stresses
# f(sigma) = F*(sigma_2 - sigma_3)^2 + G*(sigma_3 - sigma_1)^2 + H*(sigma_1 - sigma_2)^2
val = F * (s2 - s3)**2 + G * (s3 - s1)**2 + H * (s1 - s2)**2

verts, faces, normals, values = measure.marching_cubes(val, level=1.0, 
                                                       spacing=(2*lim/(res-1), 2*lim/(res-1), 2*lim/(res-1)))
# Adjust origen
verts -= lim

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

mesh = ax.plot_trisurf(verts[:, 0], verts[:, 1], faces, verts[:, 2], 
                       cmap='viridis', alpha=0.8, edgecolor='none')

line_val = np.linspace(-lim, lim, 100)
ax.plot(line_val, line_val, line_val, color='red', linestyle='--', linewidth=2, label='"Hydrostatic axis" ($\sigma_1=\sigma_2=\sigma_3$)')

# Axis setting
ax.set_xlabel(r'$\sigma_1$', fontsize=12)
ax.set_ylabel(r'$\sigma_2$', fontsize=12)
ax.set_zlabel(r'$\sigma_3$', fontsize=12)
ax.set_title('Hills Yield Surface', fontsize=14)

ax.set_xlim([-lim, lim])
ax.set_ylim([-lim, lim])
ax.set_zlim([-lim, lim])

# Perspective
ax.view_init(elev=22, azim=32, roll=1)
ax.legend()

plt.colorbar(mesh, ax=ax, shrink=0.5, aspect=10, label='Contour intensity')
plt.savefig(f"Superficie3D.png", dpi=300)
plt.close()
shutil.move(f"Superficie3D.png", os.path.join(reports_dir_cell, f"Superficie3D.png")) # move to correct directory
