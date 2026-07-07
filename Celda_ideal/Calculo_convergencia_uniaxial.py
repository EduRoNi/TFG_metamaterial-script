'''
PyMAPDL script to determine the effective mechanical properties of metamaterial RVEs via tensile testing in the three principal directions (E and nu).
'''
import matplotlib.pyplot as plt
import h5py
import numpy as np
from ansys.mapdl.core import launch_mapdl
from pathlib import Path
import os
import shutil
import sys

import openpyxl
import pandas as pd
import time

from Funciones import *

######################################################################
##                         TIME CONTROL                             ##
######################################################################
inicio_programa = time.time() # initial reference
######################################################################
##                  FILES NAMES AND DIRECTORY PATH                  ##
######################################################################

path = os.getcwd() # current working directory path

if len(sys.argv) >= 3:
    file_name = sys.argv[1]
    reports_dir = sys.argv[2]
    temporal = sys.argv[3]
reports_dir_cell = os.path.join(reports_dir, file_name)
if os.path.exists(reports_dir_cell) == False:
    os.makedirs(reports_dir_cell) # create directory to save final reports, in case it does not exist
archivo = os.path.abspath(f"{file_name}.x_t") # find the file in the current working dirextory path

tolerancia = 5e-3 # tolerance in GPa
diff = 1
structure = 4 # initialiased with 4x4x4

evo_EX = [] # EX evolution
evo_EY = [] # EY evolution
evo_EZ = [] # EZ evolution
num_celdas = [] # number of cells evolution
intento = 0 
try: 
    while diff > tolerancia:
        intento += 1
        working_dir = os.path.join(temporal, f"{file_name}.intent{intento}") # Create a temporary directory to store and delete obsolete files
        mapdl = launch_mapdl(run_location=working_dir, override=True)

        print(f"Try {intento}")
        mapdl.finish()
        mapdl.clear()
        ######################################################################
        ##                          INITIALISATION                          ##
        ######################################################################

        n = round(structure/2+0.1) # reduce number of cells for symmetry
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

        # Dimensions of complete structure (symmetry)
        dist_x = n*Lx 
        dist_y = n*Ly
        dist_z= n*Lz
        
        # Calculates effective rho
        if intento == 1:
            rho_eff = rho(mapdl, Lx, Ly, Lz, 4430)
            print(f"Effective rho = {rho_eff} kg/m^3")
        
        ##-----MESH
        mapdl.et("", "SOLID187") # SOLID187 is tetrahedron of 10 nodes
        ele_size = Lx/15 # element size using cearactheristic length
        mapdl.allsel()
        mapdl.mshkey(key="0") #0=Use free meshing (better for convergence, if not, mesh could not be created sometimes)
        mapdl.lesize(nl1="ALL", size = str(ele_size))
        mapdl.vmesh("ALL")

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
        mapdl.asel("S", "LOC", "Z", dist_z, kswp=1) # select areas normal to Z at coordenate z=dist_z
        mapdl.cm("symm_Z", "AREA") # save these areas as 'symm_Z'
        mapdl.allsel()

        mapdl.asel("S", "LOC", "Y", dist_y, kswp=1) # select area normal to Y at coordenate y=dist_y
        mapdl.cm("symm_Y", "AREA") # save these areas as 'symm_Y'
        mapdl.allsel()

        mapdl.asel("S", "LOC", "X", dist_x, kswp=1) # select area normal to X at coordenate x=dist_x
        mapdl.cm("symm_X", "AREA") # save these areas as 'symm_X'
        mapdl.allsel()
        print("Symmetry generated")

        mapdl.asel("S", "LOC", "Z", 0, kswp=1) # select areas normal to Z at coordenate z=0
        mapdl.cm("CC_Z", "AREA") # save these areas as 'CC_Z'
        A_z = dist_x*dist_y
        mapdl.allsel()
        mapdl.asel("S", "LOC", "X", 0, kswp=1) # select area normal to X at coordenate x=0
        mapdl.cm("CC_X", "AREA") # save these areas as 'CC_X'
        A_x = dist_y*dist_z
        mapdl.allsel()
        mapdl.asel("S", "LOC", "Y", 0, kswp=1) # select area normal to Y at coordenate y=0
        mapdl.cm("CC_Y", "AREA") # save these areas as 'CC_Y'
        A_y = dist_z*dist_x
        mapdl.allsel()

        # Defining symmetry
        mapdl.da("symm_X", "UX", value1=0)
        mapdl.da("symm_Y", "UY", value1=0)
        mapdl.da("symm_Z", "UZ", value1=0)

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
        mapdl.emodif("ALL", "MAT", 6)

        ##----- SAVE GEOMETRY
        # mapdl.igesout(fname=file_geometry, ext="igs", att="1")
        # print(f"Geometry file saved as: {file_geometry}")

        mapdl.allsel()
        mapdl.finish() # Close pre-processing stage

        
        # BEGIN ELASTIC ANALYSIS FOR EACH PRINCIPAL DIRECTION
        estruc_long = {"X":dist_x, "Y": dist_y, "Z": dist_z}
        celda_long = {"X":Lx, "Y": Ly, "Z": Lz}
        areas = {"X": A_x, "Y": A_y, "Z": A_z}

        Direcciones = ["X", "Y", "Z"]
        for dir in Direcciones:
            delta_L = celda_long[dir] * 0.01 # 1% cell length
            solution_1D(mapdl, f"CC_{dir}", f"U{dir}", delta_L, 0) # solver
            
            ######################################################################
            ##                 /POST1 o /POST26 |POST-PROCESS|                  ##
            ######################################################################
            mapdl.post1() # Run post-processing stage

            ##----- RESULTS
            n_sets = mapdl.post_processing.nsets
            mapdl.set(lstep="LAST", sbstep="LAST") # Analyse last results
            mapdl.allsel()

            # Deformed structure
            # mapdl.run("/RGB,INDEX,100,100,100,0")   # white background
            # mapdl.run("/RGB,INDEX,0,0,0,15")        # black text
            # mapdl.run("/REPLOT")                    
            # mapdl.plnsol("S", "EQV", 2) # deformed structure

            # Geometry
            node_ID = mapdl.mesh.nnum # nodes IDs
            node_coord = mapdl.mesh.nodes # nodes coordenates
            nodes = np.zeros((len(node_ID), 4)) # matrix = [ID x y z] for the HDF5 file
            nodes [:,0] = node_ID.T
            nodes [:,1:] = node_coord

            elem_ID = mapdl.mesh.enum # elements IDs
            elem_nodes = np.array(mapdl.mesh.elem)[:,10:20] # elements nodes
            elems = np.zeros((len(elem_ID),11)) # matrix = [ID N1 N2 N3 N4 ...] for the HDF5 file
            elems [:,0] = elem_ID.T
            elems [:,1:] = elem_nodes

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
            mapdl.cmsel("S", "symm_X", "AREA")
            mapdl.cmsel("A", "symm_Y", "AREA")
            mapdl.cmsel("A", "symm_Z", "AREA")
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

            ##------------------------------------------------------------------------------------------------##
            file_results = rf"{reports_dir_cell}\{file_name}_results_{structure}x{structure}x{structure}_{dir}"
            #----- EXPORT DATA TO HDF5 FIEL (.h5)
            with h5py.File(file_results, "w", track_order=True) as f:
                # track_order=True Save data in order of creation instead of default alphabetical order
                # Metadatos
                f.attrs["modelo"]    = file_name
                f.attrs["n_nodos"]   = len(node_ID)
                f.attrs["n_elementos"] = len(elem_ID)

                # Geometry
                geo = f.create_group("Geometry")
                geo.create_dataset("Nodes",    data=nodes) # (ID x y z)
                geo.create_dataset("Elements", data=elems) # (ID N1 N2 N3 N4 ...)
                geo.create_dataset("Constraints", data=cst) # (node_ID, "DOF", REAL_value, IMAG_value)

            for i in range(1, n_sets+1):
                mapdl.set(sbstep=str(i))
                mapdl.allsel()
                displacement,stress,strain,plastic=results(mapdl, node_ID, elem_ID)
                escritura_h5(file_results, displacement, stress, strain, plastic, i)
            ##------------------------------------------------------------------------------------------------##
            # shutil.move(file_results, os.path.join(reports_dir_cell, f"{file_name}_results_{structure}x{structure}x{structure}")) # move the file for AI to correct directory

            print(f"Save as: {file_results}")
            
            ######################################################################
            ##                             CALCULUS                             ##
            ######################################################################
            mapdl.set(lstep="LAST", sbstep="LAST") # Analyse last results
            mapdl.allsel()
            mapdl.finish() # Close post-processing stage

            displacement,stress,strain,plastic=results(mapdl, node_ID, elem_ID)
            mapdl.cmsel("S", f"CC_{dir}", "AREA")
            mapdl.nsla("S", 1)
            nodo_ref = mapdl.mesh.nnum[1] # reference node
            pos = np.where(displacement[:,0] == nodo_ref)[0]
            
            if dir == "X":
                L = dist_x
                # delta_L = abs(displacement[int(pos[0]), 1]) # actual displacement 
                E,s,e, nu_x, fff = prop_mec(mapdl, delta_L, L, A_x, dir) # Young Modulus in x direction [Pa]
                evo_EX.extend([E/1e9]) # [GPa]
            if dir == "Y":
                L = dist_y
                # delta_L = abs(displacement[int(pos[0]), 2])
                E,s,e, nu_y, fff= prop_mec(mapdl, delta_L, L, A_y, dir) # Young Modulus in y direction
                evo_EY.extend([E/1e9])
            if dir == "Z":
                L = dist_z
                # delta_L = abs(displacement[int(pos[0]), 3])
                E,s,e, nu_z, fff= prop_mec(mapdl, delta_L, L, A_z, dir) # Young Modulus in z direction
                evo_EZ.extend([E/1e9])

        # Check convergence tolerance
        if intento != 1:  
            diff = abs(E-E_viejo)/1e9 
            print(f"diff={diff}") 
        
        if diff > tolerancia:    
            structure = structure + 2
            E_viejo = E        
        
        # elif diff < tolerancia and intento < 5: 
        #     structure = structure + 2
        #     E_viejo = E
        #     diff = 1

        num_celdas.extend([n])

        # Close session and clean up temporary files
        mapdl.exit()
        shutil.rmtree(working_dir)

    print(f"Convergence at try {intento} with an structure of {structure} cells, reduced to {n} cells")

    # Determination of "n" in each principal direction
    n_exp = np.array([0.0,0.0,0.0]) # [n_exp_X, n_exp_Y, n_exp_Z]
    n_exp[0] = gibson_ashby(4430, 110, rho_eff, evo_EX[-1])
    n_exp[1] = gibson_ashby(4430, 110, rho_eff, evo_EY[-1])
    n_exp[2] = gibson_ashby(4430, 110, rho_eff, evo_EZ[-1])
    print(f"Effective rho = {rho_eff} kg/m^3")
    print(f"It has been determined that: \nX direction the structure obtain n={n_exp[0]}\nY direction the structure obtain n={n_exp[1]}\nZ direction the structure obtain n={n_exp[2]}")

except Exception as e: # if any error occurs, close all
    print(e)
    mapdl.finish()
    mapdl.exit()
    shutil.rmtree(working_dir)
    # if os.path.exists(os.path.join(path, f"{file_name}.x_t")) == True:
    #     os.remove(os.path.join(path, f"{file_name}.x_t")) # delete model from current working directory

finally:     
    fin_programa = time.time()
    # Save metadata
    meta = pd.DataFrame({
        "Parámetro": ["Model", "Final structure cells", "Convergence try",
                    "Effective density (kg/m3)",
                    "EX (GPa)", "EY (GPa)", "EZ (GPa)", "Nu_XY", "Nu_XZ", "Nu_YZ",
                    "n_exp X", "n_exp Y", "n_exp Z"],
        "Valor": [file_name, structure, intento,
                np.round(rho_eff, 2),
                evo_EX[-1], evo_EY[-1], evo_EZ[-1], nu_x[0], nu_x[1], nu_y[1], 
                n_exp[0], n_exp[1], n_exp[2]]
        })

    # Save extracted data to file
    datos_elast = {
        "Iteración": list(range(1, intento + 1)),
        "Número de celdas (n)": num_celdas,
        "Evolución de EX (GPa)": evo_EX,
        "Evolución de EY (GPa)": evo_EY,
        "Evolución de EZ (GPa)": evo_EZ
        }
    
    df_elast = pd.DataFrame(datos_elast)
    with pd.ExcelWriter(f"reporte_{file_name}.xlsx", engine='openpyxl') as writer:
        meta.to_excel(writer, sheet_name="Info", index=False)
        df_elast.to_excel(writer, sheet_name="Convergence E", index=False)

    # files control 
    shutil.move(f"reporte_{file_name}.xlsx", os.path.join(reports_dir_cell, f"reporte_{file_name}.xlsx")) # move the file to correct directory

    print("File XLSX sucessfully saved.")
    print(f"Total run time={fin_programa-inicio_programa}")
    print(f"EX = {evo_EX} \nEY = {evo_EY} \nEZ = {evo_EZ}")