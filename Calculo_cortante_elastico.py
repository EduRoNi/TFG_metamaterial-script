'''
PyMAPDL script to determine the effective mechanical properties of metamaterial RVEs via simple shear testing in the three principal directions (G).
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

archivo = os.path.abspath(f"{file_name}.x_t") # find the file
reports_dir_cell = os.path.join(reports_dir, file_name)
file_path = os.path.join(reports_dir_cell, f"reporte_{file_name}.xlsx")

# Read file
meta = pd.read_excel(file_path, sheet_name="Info")
structure = meta.loc[meta["Parámetro"] == "Final structure cells", "Valor"].values[0]

G = []

try:
    working_dir = os.path.join(temporal, f"sim_G18") # Create a temporary directory to store and delete obsolete files
    mapdl = launch_mapdl(run_location=working_dir, override=True)
    mapdl.finish()
    mapdl.clear()

    ######################################################################
    ##                          INITIALISATION                          ##
    ######################################################################
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

    # Dimensions of complete structure (symmetry)
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

    print("Geometry created")

    ##----- GENERATE SYMMETRY
    print("Generating symmetry")

    mapdl.allsel() # always select all geometry for maintain it 'active'
    mapdl.asel("S", "LOC", "Z", dist_z, kswp=1) # select areas normal to Z at coordenate z=dist_z
    mapdl.cm("disp_Z", "AREA") # save these areas as 'disp_Z'
    mapdl.allsel()

    mapdl.asel("S", "LOC", "Y", dist_y, kswp=1) # select area normal to Y at coordenate y=dist_y
    mapdl.cm("disp_Y", "AREA") # save these areas as 'disp_Y'
    mapdl.allsel()

    mapdl.asel("S", "LOC", "X", dist_x, kswp=1) # select area normal to X at coordenate x=dist_x
    mapdl.cm("disp_X", "AREA") # save these areas as 'disp_X'
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
    Direcciones = ["X", "Y", "Z"] # these directions indicates prescribed displecement
    for dir in Direcciones:
        delta_L = Lz*0.01 #1% cell length
        if dir == "X": # For displacement in X, plane normal to Z (XZ)
            solution_G(mapdl, "CC_Z", "Disp_Z", f"U{dir}", delta_L, 0)
        if dir == "Z": # For displacement in Z, plane normal to Y (YZ)
            solution_G(mapdl, "CC_Y", "Disp_Y", f"U{dir}", delta_L, 0)
        if dir == "Y": # For displacement in Y, plane normal to X (XY)
            solution_G(mapdl, "CC_X", "Disp_X", f"U{dir}", delta_L, 0)
        ######################################################################
        ##                 /POST1 o /POST26 |POST-PROCESS|                  ##
        ######################################################################
        mapdl.post1() # Run post-processing stage

        ##----- RESULTS
        n_sets = mapdl.post_processing.nsets
        mapdl.set(lstep="LAST", sbstep="LAST") # analyse last result
        mapdl.allsel()

        # Deformed structure
        # mapdl.run("/RGB,INDEX,100,100,100,0")   # white background
        # mapdl.run("/RGB,INDEX,0,0,0,15")        # black text
        # mapdl.run("/REPLOT")                    
        # mapdl.view(1, 0,1,0)
        # mapdl.plnsol("S", "EQV", 2) # deformed structure

        # Geometry
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

        file_results = rf"{reports_dir_cell}\{file_name}_results_shear_{structure}x{structure}x{structure}_{dir}"
        #----- EXPORT DATA TO HDF5 FIEL (.h5)
        with h5py.File(file_results, "w", track_order=True) as f:
            # track_order=True Save data in order of creation instead of default alphabetical order
            # Metadatos
            f.attrs["modelo"]    = file_name
            f.attrs["n_nodos"]   = len(node_ID)
            f.attrs["n_elementos"] = len(elem_ID)

            # Geometría
            geo = f.create_group("Geometry")
            geo.create_dataset("Nodes",    data=nodes) # (ID x y z)
            geo.create_dataset("Elements", data=elems) # (ID N1 N2 N3 N4 ...)
            geo.create_dataset("Constraints", data=cst) 

        for i in range(1, n_sets+1):
            mapdl.set(sbstep=str(i))
            mapdl.allsel()
            displacement,stress,strain,plastic=results(mapdl, node_ID, elem_ID)
            escritura_h5(file_results, displacement, stress, strain, plastic, i)


        if dir == "X": # G_XZ
            g, uat, ammag, FF = prop_mec_G(mapdl, delta_L, dist_z, dist_x*dist_y, "Z", dir)
        if dir == "Z": # G_YZ
            g, uat, ammag, FF = prop_mec_G(mapdl, delta_L, dist_y, dist_x*dist_z, "Y", dir)
        if dir == "Y": # G_XY
            g, uat, ammag, FF = prop_mec_G(mapdl, delta_L, dist_x, dist_y*dist_z, "X", dir)
        
        G.extend([g/1e9]) # [Pa]    

    print(f"For {file_name} the next G values are obtained:")
    print(f"Gxy = {G[1]} GPa \nGxz = {G[0]} GPa \nGyz = {G[2]} GPa")

    # File control
    meta_new = pd.concat([meta, pd.DataFrame({
        "Parámetro": ["Gxy (GPa)", "Gxz (GPa)", "Gyz (GPa)"],
        "Valor": [G[1], G[0], G[2]]
    })], ignore_index=True)

    # Overwrite to include new information
    with pd.ExcelWriter(file_path, engine='openpyxl', mode='a', if_sheet_exists='replace') as writer:
        meta_new.to_excel(writer, sheet_name="Info", index=False)

finally:
    #close and clear all
    mapdl.finish()
    mapdl.exit()
    shutil.rmtree(working_dir)

