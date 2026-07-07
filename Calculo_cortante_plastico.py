'''
PyMAPDL script to determine the elastoplastic behavior of metamaterial RVEs via simple shear testing in the three principal directions.
'''
import matplotlib.pyplot as plt
import h5py
import numpy as np
from ansys.mapdl.core import launch_mapdl
from pathlib import Path 
import os
import shutil
import sys

import pandas as pd
import time

from Funciones import *

######################################################################
##                  FILES NAMES AND DIRECTORY PATH                  ##
######################################################################
path = os.getcwd() # current working directory path
if len(sys.argv) >= 3:
    file_name = sys.argv[1]
    reports_dir = sys.argv[2]
    temporal = sys.argv[3]

archivo = os.path.abspath(f"{file_name}.x_t") # find the file
reports_dir_cell = os.path.join(reports_dir, file_name)
file_path = os.path.join(reports_dir_cell, f"reporte_{file_name}.xlsx")

inicio_programa = time.time()

# Read file to extract infromation about the cell
meta = pd.read_excel(file_path, sheet_name="Info")
structure = meta.loc[meta["Parámetro"] == "Final structure cells", "Valor"].values[0] #complete structure
Gxy = meta.loc[meta["Parámetro"] == "Gxy (GPa)", "Valor"].values[0]
Gxz = meta.loc[meta["Parámetro"] == "Gxz (GPa)", "Valor"].values[0]
Gyz = meta.loc[meta["Parámetro"] == "Gyz (GPa)", "Valor"].values[0]
G_array = [Gxy*1e9, Gxz*1e9, Gyz*1e9] # G_array = [GXY, GXZ, GYZ] [Pa]
print(G_array)
direcciones = ["X", "Y", "Z"]
try:
    for j in range(0,len(direcciones)):
        try:
            dir = direcciones[j]
            print(f"ANALYSE IN {dir}")
            working_dir = os.path.join(temporal, f"sim_G_plast_{dir}") # Create a temporary directory to store and delete obsolete files
            mapdl = launch_mapdl(run_location=working_dir, override=True)
            mapdl.finish()
            mapdl.clear()

            ######################################################################
            ##                          INITIALISATION                          ##
            ######################################################################
            file_results = rf"{reports_dir_cell}\{file_name}_plastic_results_shear_{structure}x{structure}x{structure}_{dir}"
            n = structure
            ##----- READ FILE
            mapdl.parain(name=file_name, extension="x_t", path=archivo, fmt="0", scale="0")
            print(f"Model {file_name} has been opened")

            ######################################################################
            ##                       /PREP7 |PRE-PROCESS|                       ##
            ######################################################################
            mapdl.prep7() # Run pre-processing stage
            mapdl.units("SI") # International System of Units

            ##----- GEOMETRY DATA
            keypoints = mapdl.geometry.get_keypoints(return_as_array=True) #Return the keypoints coordinates as a numpy array = [[x, y, z], [x, y, z], ...]
            x_max = max(keypoints[:,0]) # max x coord
            x_min = min(keypoints[:,0]) # min x coord

            y_max = max(keypoints[:,1]) # max y coord
            y_min = min(keypoints[:,1]) # min y coord

            z_max = max(keypoints[:,2]) # max z coord
            z_min = min(keypoints[:,2]) # min z coord

            # Dimensions of a cell
            Lx = x_max - x_min
            Ly = y_max - y_min
            Lz = z_max - z_min

            # Dimensions of complete structure
            dist_x = n*Lx 
            dist_y = n*Ly
            dist_z= n*Lz

            ###############################
            ##-----MESH
            mapdl.et("", "SOLID187") # SOLID187 is tetrahedron of 10 nodes
            ele_size = Lx/15 # element size using cearactheristic length
            mapdl.allsel()
            mapdl.mshkey(key="0") #0=Use free meshing (better for convergence, if not, mesh could not be created sometimes)
            # mapdl.esize(size=str(ele_size)) 
            # mapdl.smrtsize(1)
            mapdl.lesize(nl1="ALL", size = str(ele_size))
            mapdl.vmesh("ALL")
            ##############################

            ##----- GENERATE NEW STRUCTURE
            print("Geometry creation")
            mapdl.allsel()
            mapdl.vgen(itime=str(n), nv1="ALL", dz=str(Lz), noelem="0", imove="0") # line of cells in z direction

            mapdl.allsel()
            mapdl.vgen(itime=str(n), nv1="ALL", dy=str(Ly), noelem="0", imove="0") # plane of cells in zy

            mapdl.allsel()
            mapdl.vgen(itime=str(n), nv1="ALL", dx=str(Lx), noelem="0", imove="0") # complete geometry
            mapdl.allsel()

            ##################################
            # Merge coincident nodes to connect adjacent geometries
            mapdl.nummrg("NODE", ele_size * 0.01)
            mapdl.allsel()
            ##################################

            ##----- GENERATE SYMMETRY
            print("Generating symmetry")
            mapdl.allsel() # always select all geometry for maintain it 'active'
            mapdl.asel("S", "LOC", "Z", dist_z) # select areas normal to Z at coordenate z=dist_z
            mapdl.cm("disp_Z", "AREA") # save these areas as 'disp_Z'
            mapdl.allsel()

            mapdl.asel("S", "LOC", "Y", dist_y) # select area normal to Y at coordenate y=dist_y
            mapdl.cm("disp_Y", "AREA") # save these areas as 'disp_Y'
            mapdl.allsel()

            mapdl.asel("S", "LOC", "X", dist_x) # select area normal to X at coordenate x=dist_x
            mapdl.cm("disp_X", "AREA") # save these areas as 'disp_X'
            mapdl.allsel()
            print("Symmetry generated")

            mapdl.asel("S", "LOC", "Z", 0) # select areas normal to Z at coordenate z=0
            mapdl.cm("CC_Z", "AREA") # save these areas as 'CC_Z'
            A_z = dist_x*dist_y
            mapdl.allsel()
            mapdl.asel("S", "LOC", "X", 0) # select area normal to X at coordenate x=0
            mapdl.cm("CC_X", "AREA") # save these areas as 'CC_X'
            A_x = dist_y*dist_z
            mapdl.allsel()
            mapdl.asel("S", "LOC", "Y", 0) # select area normal to Y at coordenate y=0
            mapdl.cm("CC_Y", "AREA") # save these areas as 'CC_Y'
            A_y = dist_z*dist_x
            mapdl.allsel()

            ##----- MATERIAL DEFINITION
            # Defining Structural Steel Lineal Isotropic material as mat=1
            mapdl.mp(lab="EX", mat="1", c0=200e9) #Young Modulus [Pa]
            mapdl.mp(lab="PRXY", mat="1", c0=0.3) #Poisson ratio
            mapdl.mp(lab="DENS", mat="1", c0=7850) #Density [kg/m3]
            su = 410e6 # ultimate tensile stress [Pa]

            # Defining Structural Steel Bilineal Isotropic Hardening material as mat=2
            mapdl.mp(lab="EX", mat="2", c0=200e9) #Young Modulus [Pa]
            mapdl.mp(lab="PRXY", mat="2", c0=0.3) #Poisson ratio
            mapdl.mp(lab="DENS", mat="2", c0=7850) #Density [kg/m3]
            mapdl.tb(lab="BISO", mat="2") #BISO: Bilinear isotropic hardening using von Mises or Hill plasticity.
                                        #PLASTIC: Nonlinear plasticity (Multilinear isotropic hardening)
                                        #NLISO: Voce isotropic hardening law (or power law) for modeling nonlinear isotropic hardening using von Mises or Hill plasticity.
            mapdl.tbdata(stloc="1", c1=280e6, c2=500e6) #.tbdata(1, c1=Yield stress [Pa], c2=Plastic Tangent Modulus[Pa])
            # mapdl.tbplot("BISO", mat="2")

            #Defining Structutal Steel Power Law Nonlinear Isotropic Hardening as mat=3
            mapdl.mp(lab="EX", mat="3", c0=200e9) #Young Modulus [Pa]
            mapdl.mp(lab="PRXY", mat="3", c0=0.3) #Poisson ratio
            mapdl.mp(lab="DENS", mat="3", c0=7850) #Density [kg/m3]
            mapdl.tb("NLISO", mat=3, funcname="POWER")
            mapdl.tbdata(stloc="1", c1=250e6, c2=0.15) #.tbdata(1, c1=Initial yield stress [Pa], c2=Exponent [-])

            #Defining Structural Steel Voce Law Nonlinear Isotropic Hardening as mat=4
            mapdl.mp(lab="EX", mat="4", c0=200e9) #Young Modulus [Pa]
            mapdl.mp(lab="PRXY", mat="4", c0=0.3) #Poisson ratio
            mapdl.mp(lab="DENS", mat="4", c0=7850) #Density [kg/m3]
            mapdl.tb("NLISO", mat=4, funcname="VOCE")
            mapdl.tbdata(stloc="1", c1=250e6, c2=0, c3=450e6, c4=15) #.tbdata(1, c1=Initial yield stress [Pa], c2=Linear coefficient [Pa], c3=Exponential coefficient [Pa], c4=Exponential Saturation Parameter [-])

            #Defining Ti-6Al-4V Multilinear Isotropic Hardening as mat=5
            mapdl.mp(lab="EX", mat="5", c0=110*1e9) #Young Modulus [Pa]
            mapdl.mp(lab="PRXY", mat="5", c0=0.33) #Poisson ratio
            mapdl.mp(lab="DENS", mat="5", c0=4430) #Density [kg/m3]
            mapdl.tb("PLASTIC", mat=5, funcname="MISO")
            # Define plastic stress-strain curve properties (MISO)
            mapdl.tbpt(oper="DEFI", x1=0.000000, x2=952.490 * 1e6)
            mapdl.tbpt(oper="DEFI", x1=0.006400, x2=976.800 * 1e6)
            mapdl.tbpt(oper="DEFI", x1=0.011377, x2=994.662 * 1e6)
            mapdl.tbpt(oper="DEFI", x1=0.017777, x2=1010.522 * 1e6)
            mapdl.tbpt(oper="DEFI", x1=0.025598, x2=1027.591 * 1e6)
            mapdl.tbpt(oper="DEFI", x1=0.034842, x2=1046.823 * 1e6)
            mapdl.tbpt(oper="DEFI", x1=0.045508, x2=1066.043 * 1e6)
            mapdl.tbpt(oper="DEFI", x1=0.057596, x2=1081.125 * 1e6)
            mapdl.tbeo("NEGSLOPE", 1) # allows negative tangent slope, just in case 
            su = 1081.125 * 1e6 # ultimate tensile stress [Pa]
            plastic_lim = 0.057596 # ultimate plastic strain (failure criterion)
            
            #Defining Ti-6Al-4V Lineal Isotropic as mat=6
            mapdl.mp(lab="EX", mat="6", c0=110*1e9) #Young Modulus [Pa]
            mapdl.mp(lab="PRXY", mat="6", c0=0.33) #Poisson ratio
            mapdl.mp(lab="DENS", mat="6", c0=4430) #Density [kg/m3]

            # Element material attribute
            mapdl.allsel()
            mapdl.esel("ALL")
            mapdl.emodif("ALL", "MAT", 5)

            ##----- SAVE GEOMETRY
            # mapdl.igesout(fname=file_geometry, ext="igs", att="1")
            # print(f"Geometry file saved as: {file_geometry}")

            mapdl.allsel()
            mapdl.finish() # Close pre-processing stage

            ######################################################################
            ##            /SOLU |SOLVER| y /POST1 |POSTPROCESS|                 ##
            ######################################################################
            mapdl.filname() # Changes the Jobname for the analysis to "work_dir"
            mapdl.slashsolu() 

            if dir == "X":  # For displacement in X, plane normal to Z (XZ)
                A = A_z # [m^2]
                dist = dist_z # characteristic length
                dir2 = "Z" # Direction perpendicular to the area where the displacement is imposed
                #BOUNDARY CONDITIONS
                # FIXED
                mapdl.cmsel("S", "CC_Z", "AREA")
                mapdl.nsla("S", 1)
                mapdl.d(node="ALL", lab="ALL", value=0)

                # CONSTRAINTS LATERAL FACES
                mapdl.cmsel("S", "Disp_X", "AREA")
                mapdl.cmsel("A", "CC_X", "AREA")
                mapdl.cmsel("A", "Disp_Y", "AREA")
                mapdl.cmsel("A", "CC_Y", "AREA")
                mapdl.nsla("S", 1)
                mapdl.nsla("S", 1)
                mapdl.d(node="ALL", lab="UY", lab2="UZ", value=0)

                # CONSTRAINTS UPPER FACES
                mapdl.allsel()
                mapdl.cmsel("S", "Disp_Z", "AREA")
                mapdl.nsla("S", 1)
                mapdl.d(node="ALL", lab="UY", lab2="UZ", value=0)
            if dir == "Y": # For displacement in Y, plane normal to X (XY) 
                A = A_x
                dist = dist_x 
                dir2 = "X"
                #BOUNDARY CONDITIONS
                # FIXED
                mapdl.cmsel("S", "CC_X", "AREA")
                mapdl.nsla("S", 1)
                mapdl.d(node="ALL", lab="ALL", value=0)

                # CONSTRAINTS LATERAL FACES
                mapdl.cmsel("S", "Disp_Y", "AREA")
                mapdl.cmsel("A", "CC_Y", "AREA")
                mapdl.cmsel("A", "Disp_Z", "AREA")
                mapdl.cmsel("A", "CC_Z", "AREA")
                mapdl.nsla("S", 1)
                mapdl.d(node="ALL", lab="UX", lab2="UZ", value=0)

                # CONSTRAINTS UPPER FACE
                mapdl.allsel()
                mapdl.cmsel("S", "Disp_X", "AREA")
                mapdl.nsla("S", 1)
                mapdl.d(node="ALL", lab="UX", lab2="UZ", value=0)
            if dir == "Z": # For displacement in Z, plane normal to Y (YZ)
                A = A_y
                dist = dist_y
                dir2 = "Y"
                #BOUNDARY CONDITIONS
                # FIXED
                mapdl.cmsel("S", "CC_Y", "AREA")
                mapdl.nsla("S", 1) 
                mapdl.d(node="ALL", lab="ALL", value=0)

                # CONSTRAINTS LATERAL FACES
                mapdl.cmsel("S", "Disp_Z", "AREA")
                mapdl.cmsel("A", "CC_Z", "AREA")
                mapdl.cmsel("A", "Disp_X", "AREA")
                mapdl.cmsel("A", "CC_X", "AREA")
                mapdl.nsla("S", 1)
                mapdl.d(node="ALL", lab="UY", lab2="UX", value=0)

                # CONSTRAINTS UPPER FACE
                mapdl.allsel()
                mapdl.cmsel("S", "Disp_Y", "AREA")
                mapdl.nsla("S", 1)
                mapdl.d(node="ALL", lab="UY", lab2="UX", value=0)

            rotura = False  # Define solver stopping variable
            increment = dist*0.05 # 5% structure's characteristic length
            delta_L = 0

            mapdl.nlgeom("ON")
            mapdl.eqslv("SPARSE") 
            mapdl.nropt("FULL") #use FULL Newton-Rapshon method

            # Arrays for saving data
            tau_global = [] # global tau values
            gamma_global = [] # global gamma values
            FR = [] # reaction force
            disp_global = [] # global displacement
            load_step = 0 # load steps done
            old_substp = 0 # total substeps

            primera_vez = True

            while  rotura == False:
                mapdl.allsel()
                delta_L += increment # increment displacement
                mapdl.cmsel("S", f"Disp_{dir2}", "AREA")
                mapdl.nsla("S", 1)
                mapdl.d(node="ALL", lab=f"U{dir}", value=+delta_L)

                if primera_vez == True: # First iteration, initilise simulation
                    print("new")
                    mapdl.antype("STATIC", "NEW")
                    mapdl.rescontrol("", "ALL", -1) # Save file for resatr secuence
                    primera_vez = False
                else: # Next iteration, restart from the previous load step
                    print("restart")
                    print(load_step, new_substp)
                    mapdl.antype("", "REST", load_step, new_substp) # start from the last substep of the last step
                    mapdl.cmsel("S", f"Disp_{dir2}", "AREA")
                    mapdl.nsla("S", 1)
                    mapdl.d(node="ALL", lab=f"U{dir}", value=+delta_L)

                mapdl.nsubst(10,20,5) # substeps definition
                mapdl.autots("ON")  # automatic time stepping ON  
                mapdl.outres("ALL", "ALL") # save all substeps

                ##----- SOLVER
                print("Solver initialised")
                mapdl.allsel()
                mapdl.solve()
                mapdl.finish()

                mapdl.post1()
                n_sets = mapdl.post_processing.nsets
                mapdl.allsel()

                node_ID = mapdl.mesh.nnum # nodes IDs
                node_coord = mapdl.mesh.nodes # nodes coordenates
                nodes = np.zeros((len(node_ID), 4)) # matriz = [ID x y z] for the HDF5 file
                nodes [:,0] = node_ID.T
                nodes [:,1:] = node_coord

                elem_ID = mapdl.mesh.enum # elements IDs
                elem_nodes = np.array(mapdl.mesh.elem)[:,10:20] # elements nodes
                elems = np.zeros((len(elem_ID),11)) # matriz = [ID N1 N2 N3 N4 ...] for the HDF5 file
                elems [:,0] = elem_ID.T
                elems [:,1:] = elem_nodes

                mapdl.cmsel("S", f"Disp_{dir2}", "AREA")
                mapdl.nsla("S", 1)
                nodo_ref = mapdl.mesh.nnum[1] # reference node
                load_step += 1
                new_substp = n_sets - old_substp 
                old_substp = n_sets
                for i in range(1, new_substp+1):
                    mapdl.set(lstep= load_step, sbstep=str(i))
                    mapdl.allsel()
                    displacement,stress,strain,plastic=results(mapdl, node_ID, elem_ID)
                    escritura_h5(file_results, displacement, stress, strain, plastic, str(load_step)+str(i))
                    sigma_local = max(stress[:,7])
                    plastic_local = max(plastic[:,1])
                    ##########--------------------------################----------------------###########
                    eppl = plastic[:,1] # equivalent plastic strain
                    eppl_acc = plastic[:,2] # accumulated plastic strain
                    plastic_values = eppl[eppl > plastic_lim]
                    plastic_values_acc = eppl_acc[eppl_acc > plastic_lim] 
                    # print(plastic_values)
                    # print(plastic_values_acc)
                    ##########--------------------------################----------------------###########
                    pos = np.where(displacement[:,0] == nodo_ref)[0]
                    delta = abs(displacement[int(pos[0]), j+1]) # actual displacement
                    
                    GG, t, g, f= prop_mec_G(mapdl, delta, dist, A, dir2, dir)
                    print(f"For a total displacement of {delta} m, a local sigma = {sigma_local} is achieved. Shear modulus = {GG/1e9} GPa")
                    tau_global.extend([t])
                    gamma_global.extend([g])
                    FR.extend([f])
                    disp_global.extend([delta])
                    if len(plastic_values_acc) > len(elem_ID) * 0.05:
                        print(f"The plastic strain limit has been exceeded (plastic limit = {plastic_lim}) in {len(plastic_values)} elements exceeding the 5% of elements damage with a global displacement of {delta} m (theoretically {delta_L} m).")
                        break
                if len(plastic_values_acc) > len(elem_ID) * 0.05:
                    break

                mapdl.finish()
                mapdl.slashsolu() # return to SOLVER for new iteration

            # Restrictions and boundary conditions
            dof_labels = ["UX", "UY", "UZ"]
            records = []

            # Restricted DOF
            constrained_sets = {}
            for dof in dof_labels:
                mapdl.allsel()
                mapdl.nsel("S", "D", dof)
                constrained_sets[dof] = set(mapdl.mesh.nnum)

            # Matrix creation
            mapdl.allsel()
            mapdl.cmsel("S", "disp_X", "AREA")
            mapdl.cmsel("A", "disp_Y", "AREA")
            mapdl.cmsel("A", "disp_Z", "AREA")
            mapdl.cmsel("A", "CC_X", "AREA")
            mapdl.cmsel("A", "CC_Y", "AREA")
            mapdl.cmsel("A", "CC_Z", "AREA")
            mapdl.nsla("S", 1) # select only the nodes from the boundaries
            for n in mapdl.mesh.nnum:
                row = [n]
                for dof in dof_labels:
                    if n in constrained_sets[dof]:
                        val = mapdl.get_value(entity="NODE", entnum=n, item1="D", it1num=dof)
                    else:
                        val = 1e30  # Free DOF
                    row.append(val)
                records.append(row)

            cst = np.array(records) 

            #----- EXPORT DATA TO HDF5 FIEL (.h5)
            with h5py.File(file_results, "a", track_order=False) as f:
                # track_order=False Save data in alphabetical order
                # Metadatos
                f.attrs["modelo"]    = file_name
                f.attrs["n_nodos"]   = len(node_ID)
                f.attrs["n_elementos"] = len(elem_ID)

                # Geometry
                geo = f.create_group("Geometry")
                geo.create_dataset("Nodes",    data=nodes) # (ID x y z)
                geo.create_dataset("Elements", data=elems) # (ID N1 N2 N3 N4 ...)
                geo.create_dataset("Constraints", data=cst) # (node_ID, Ux, Uy, Uz)

            mapdl.finish()

            tau_global = np.array(tau_global)
            gamma_global = np.array(gamma_global) 
            FR = np.array(FR)
            disp_global = np.array(disp_global)

            G = G_array[j] # [Pa]
            S, inte, Miso = plast(G, tau_global, gamma_global)
            print(f"Shear modulus= {G/1e9} GPa, Yield stress S_{dir}{dir2} = {S/1e6} MPa" )
            print(f"Max internal stress= {max(stress[:,7])/1e6} MPa")

            # Graphs
            # Settings
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

            # stress-strain curve
            fig, ax = plt.subplots(figsize=(7, 6))

            ax.plot(gamma_global, tau_global/1e6, color='red', linewidth=2, label='Stress-Strain') 
            
            ax.scatter(gamma_global[inte], S/1e6, marker='o', facecolors='none', edgecolors='black', 
                    s=60, linewidths=1.2, zorder=5,
                    label=f"Yield stress {np.round(S/1e6, 2)} MPa") 
            
            ax.scatter(gamma_global[-1], tau_global[-1]/1e6, marker='x', color='black', 
                    s=60, linewidths=1.5, zorder=5,
                    label=f"Failure {np.round(tau_global[-1]/1e6, 2)} MPa") 
            
            # Labels and titles
            ax.set_title(f'Shear stress-strain {dir}{dir2}', fontsize=14, pad=10)
            ax.set_xlabel('Strain ($\gamma$)', fontsize=13, labelpad=8)
            ax.set_ylabel('\u03c4 [MPa]', fontsize=13, labelpad=8)
            
            ax.tick_params(top=True, right=True, direction='in')
            
            # Legend
            ax.legend(loc='best', edgecolor='black', framealpha=1)
            
            ax.grid(True, linestyle='-', linewidth=0.5, color='lightgray')
        
            plt.tight_layout()
            
            # Save and close
            plt.savefig(f"Curva_TTangencialDAngular_{dir}{dir2}.png", dpi=300)
            plt.close()
            shutil.move(f"Curva_TTangencialDAngular_{dir}{dir2}.png", os.path.join(reports_dir_cell, f"Curva_TTangencialDAngular_{dir}{dir2}.png")) # move to correct directory


            # force-displacement curve
            fig, ax = plt.subplots(figsize=(7, 6))

            ax.plot(disp_global, FR/1e3, color='red', linewidth=2, label='Force-Displacement') 
            
            ax.scatter(disp_global[-1], FR[-1]/1e3, marker='x', color='black', 
                    s=60, linewidths=1.5, zorder=5,
                    label=f"Failure {np.round(FR[-1]/1e3, 2)} kN") 
            
            # Labels and titles
            ax.set_title(f'Force-Displacement {dir}{dir2}', fontsize=14, pad=10)
            ax.set_xlabel(f'Displacement {dir} [mm]', fontsize=13, labelpad=8)
            ax.set_ylabel(f'Force {dir} [kN]', fontsize=13, labelpad=8)
            
            ax.tick_params(top=True, right=True, direction='in')
            
            # Legend
            ax.legend(loc='best', edgecolor='black', framealpha=1)
            
            ax.grid(True, linestyle='-', linewidth=0.5, color='lightgray')
        
            plt.tight_layout()
            
            # Save and close
            plt.savefig(f"Curva_FuerzaDesplaz_{dir}{dir2}.png", dpi=300)
            plt.close()
            shutil.move(f"Curva_FuerzaDesplaz_{dir}{dir2}.png", os.path.join(reports_dir_cell, f"Curva_FuerzaDesplaz_{dir}{dir2}.png")) # move to correct directory

            # File control
            # Meta
            meta = pd.read_excel(file_path, sheet_name="Info")
            meta_new = pd.concat([meta, pd.DataFrame({
                "Parámetro": [f"Límite elast_{dir}{dir2} (MPa)", f"Tensión tangencial de rotura_{dir}{dir2} (MPa)", f"Tensión interna máx (MPa)"],
                "Valor": [S/1e6, tau_global[-1]/1e6, sigma_local/1e6]
            })], ignore_index=True)

            # Overwrite to include new information
            with pd.ExcelWriter(file_path, engine='openpyxl', mode='a', if_sheet_exists='replace') as writer:
                meta_new.to_excel(writer, sheet_name="Info", index=False)

                # Data
                datos_plast = {
                    f"Tau_{dir}{dir2} (Pa)": tau_global,
                    f"Gamma_{dir}{dir2}": gamma_global,
                    f"Displacement_{dir} (mm)": disp_global,
                    f"Force_{dir} (N)": FR
                }
                df_plast = pd.DataFrame(datos_plast)
                df_plast.to_excel(writer, sheet_name=f"Plasticidad_{dir}{dir2}", index=False)

            print("File XLSX sucessfully saved.")

            print(f"FINISHED ANALYSIS IN {dir}{dir2}")
            mapdl.finish()
            mapdl.exit()
            shutil.rmtree(working_dir)
        
        except Exception as e:
            mapdl.finish()
            mapdl.exit()
            shutil.rmtree(working_dir)

            print(f"Convergence error during load step {load_step}. Considered as structural collapse due to buckling.")
            print(f"Internal critical load {max(stress[:,7])/1e6} MPa, with a global displacement of {delta} m (theoretically {delta_L} m)")

            tau_global = np.array(tau_global)
            gamma_global = np.array(gamma_global) 

            G = G_array[j] # [Pa]
            S, inte, Miso = plast(G, tau_global, gamma_global)
            print(f"Shear modulus= {G/1e9} GPa, Yield stress S_{dir}{dir2} = {S/1e6} MPa" )
            print(f"Max internal stress= {max(stress[:,7])/1e6} MPa")

            # Graphs
            # Settings
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

            # stress-strain curve
            fig, ax = plt.subplots(figsize=(7, 6))

            ax.plot(gamma_global, tau_global/1e6, color='red', linewidth=2, label='Stress-Strain') 
            
            ax.scatter(gamma_global[inte], S/1e6, marker='o', facecolors='none', edgecolors='black', 
                    s=60, linewidths=1.2, zorder=5,
                    label=f"Yield stress {np.round(S/1e6, 2)} MPa") 
            
            ax.scatter(gamma_global[-1], tau_global[-1]/1e6, marker='x', color='black', 
                    s=60, linewidths=1.5, zorder=5,
                    label=f"Failure {np.round(tau_global[-1]/1e6, 2)} MPa") 
            
            # Labels and titles
            ax.set_title(f'Shear stress-strain {dir}{dir2}', fontsize=14, pad=10)
            ax.set_xlabel('Strain ($\gamma$)', fontsize=13, labelpad=8)
            ax.set_ylabel('\u03c4 [MPa]', fontsize=13, labelpad=8)
            
            ax.tick_params(top=True, right=True, direction='in')
            
            # Legend
            ax.legend(loc='best', edgecolor='black', framealpha=1)
            
            ax.grid(True, linestyle='-', linewidth=0.5, color='lightgray')
        
            plt.tight_layout()
            
            # Save and close
            plt.savefig(f"Curva_TTangencialDAngular_{dir}{dir2}.png", dpi=300)
            plt.close()
            shutil.move(f"Curva_TTangencialDAngular_{dir}{dir2}.png", os.path.join(reports_dir_cell, f"Curva_TTangencialDAngular_{dir}{dir2}.png")) # move to correct directory


            # force-displacement curve
            fig, ax = plt.subplots(figsize=(7, 6))

            ax.plot(disp_global, FR/1e3, color='red', linewidth=2, label='Force-Displacement') 
            
            ax.scatter(disp_global[-1], FR[-1]/1e3, marker='x', color='black', 
                    s=60, linewidths=1.5, zorder=5,
                    label=f"Failure {np.round(FR[-1]/1e3, 2)} kN") 
            
            # Labels and titles
            ax.set_title(f'Force-Displacement {dir}{dir2}', fontsize=14, pad=10)
            ax.set_xlabel(f'Displacement {dir} [mm]', fontsize=13, labelpad=8)
            ax.set_ylabel(f'Force {dir} [kN]', fontsize=13, labelpad=8)
            
            ax.tick_params(top=True, right=True, direction='in')
            
            # Legend
            ax.legend(loc='best', edgecolor='black', framealpha=1)
            
            ax.grid(True, linestyle='-', linewidth=0.5, color='lightgray')
        
            plt.tight_layout()
            
            # Save and close
            plt.savefig(f"Curva_FuerzaDesplaz_{dir}{dir2}.png", dpi=300)
            plt.close()
            shutil.move(f"Curva_FuerzaDesplaz_{dir}{dir2}.png", os.path.join(reports_dir_cell, f"Curva_FuerzaDesplaz_{dir}{dir2}.png")) # move to correct directory

            # File control
            # Meta
            meta = pd.read_excel(file_path, sheet_name="Info")
            meta_new = pd.concat([meta, pd.DataFrame({
                "Parámetro": [f"Límite elast_{dir}{dir2} (MPa)", f"Tensión tangencial de rotura_{dir}{dir2} (MPa)", f"Tensión interna máx (MPa)"],
                "Valor": [S/1e6, tau_global[-1]/1e6, sigma_local/1e6]
            })], ignore_index=True)

            # Overwrite to include new information
            with pd.ExcelWriter(file_path, engine='openpyxl', mode='a', if_sheet_exists='replace') as writer:
                meta_new.to_excel(writer, sheet_name="Info", index=False)

                # Data
                datos_plast = {
                    f"Tau_{dir}{dir2} (Pa)": tau_global,
                    f"Gamma_{dir}{dir2}": gamma_global,
                    f"Displacement_{dir} (mm)": disp_global,
                    f"Force_{dir} (N)": FR
                }
                df_plast = pd.DataFrame(datos_plast)
                df_plast.to_excel(writer, sheet_name=f"Plasticidad_{dir}{dir2}", index=False)

            print("File XLSX sucessfully saved.")

            print(f"FINISHED ANALYSIS IN {dir}{dir2}")
            mapdl.finish()
            mapdl.exit()
            shutil.rmtree(working_dir)
finally:
    fin_programa = time.time()
    print(f"Total time={fin_programa-inicio_programa}")