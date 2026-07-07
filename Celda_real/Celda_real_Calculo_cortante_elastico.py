'''
PyMAPDL script to determine the effective mechanical properties of metamaterial RVEs via simple shear testing in the three principal directions (G).
Adapted to CeldillaV1 real.
'''
'''
Codigo que está modificado para poder estudiar la celda de metamaterial real Los cambios que contiene el código son:
    · A parte de leer el archivo con la geometría para sacar los datos geométricos, se le debe entregar un archivo en 
      el que el modelo esté previamente mallado, ya que el mallado automático aquí no funciona por la complejidad de la geometría
    · Como ahora lo que tengo son nodos y elementos, para hacer la estructura utilizo egen() en vez de vgen()
    · No se identifican AREAS así que directamente debo seleccionar nodos/elementos:
        -> cambiando todos los asel() por nsel()
        -> cambiando dentro de cm("AREA") por cm("NODE")
        -> cambiando dentro de lsa distintas funciones definidas (def) el seleccionar directamente nodos y no areas y luego nodos
        -> cond.cont. de simetría aplicadas directamente sobre nodos (mapdl.d()) y no sobre áreas (mapdl.da())
    · NO es posible fusionar los nodos pq las mallas no son las mismas así que es necesario usar mapdl.cpintf("ALL") que realciona 
      los GDL de los nodos cercanos sin eliminarlos (coupling). El problema surge de que no es posible aplicar luego cond.cont. a aquellos
      nodos 'coupled', así que debo aplicar mapdl.cpintf() a aqeullos nodos que no vayan a tener ninguna otra restricción (he creado una nueva def coupling() para
      usarla dentro de solution_1D_modif)
    · Para crear el archivo ya mallaod (es el ds.dat convertido a ds.cdb) era necesario hacer una simulación en ANSYS. Este archivo
      guarda las condiciones decontrono impuestas por lo que será preciso borrarlas todas para el cálculo
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
    Direcciones = ["X", "Y", "Z"]
    for dir in Direcciones:
        working_dir = os.path.join(temporal, f"CeldaReal.G_{dir}") # Create a temporary directory to store and delete obsolete files
        mapdl = launch_mapdl(run_location=working_dir, override=True)
        
        mapdl.finish()
        mapdl.clear()

        ######################################################################
        ##                          INITIALISATION                          ##
        ######################################################################
        file_results = rf"{path}\{file_name}_results_{structure}x{structure}x{structure}"
        file_geometry = rf"{path}\{file_name}_geometry_{structure}x{structure}x{structure}"
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

        # Delete and clean everything as we already have all relevant geometry data. Now open model previously meshed
        mapdl.finish()
        mapdl.clear()
        mapdl.cdread('db', 'ds.cdb')

        mapdl.slashsolu() # open solver
        mapdl.allsel()
        mapdl.ddele("ALL", "ALL") # delete all boundary conditions
        mapdl.finish()

        mapdl.prep7() # open pre-processor
        ##----- GENERATE NEW STRUCTURE USING ELEMENTS (NO SOLIDS/VOLUMES)
        print("Geometry creation")
        mapdl.allsel()
        mapdl.egen(itime=str(n), ninc= mapdl.mesh.n_node, iel1="ALL", dz=Lz)
        # mapdl.vgen(itime=str(n), nv1="ALL", dz=str(Lz), noelem="0", imove="0") # line of cells in z direction

        mapdl.allsel()
        mapdl.egen(itime=str(n), ninc= mapdl.mesh.n_node, iel1="ALL", dy=Ly)
        # mapdl.vgen(itime=str(n), nv1="ALL", dy=str(Ly), noelem="0", imove="0") # plane of cells in zy

        mapdl.allsel()
        mapdl.egen(itime=str(n), ninc= mapdl.mesh.n_node, iel1="ALL", dx=Lx)
        # mapdl.vgen(itime=str(n), nv1="ALL", dx=str(Lx), noelem="0", imove="0") # complete geometry
        mapdl.allsel()
        print("Geometry created")

        # ----- MERGE COINCIDENTS NODES
        ele_size = Lx/15
        for i in range(1,n):
            # X normal plane interface
            tol = 1e-6
            interface_x = i*Lx
            mapdl.nsel("S", "LOC", "X", interface_x) # select
            mapdl.nummrg(label="NODE", toler=tol) # merge
            mapdl.allsel()

            # Y normal plane interface
            interface_y = i*Ly
            mapdl.nsel("S", "LOC", "Y", interface_y)
            mapdl.nummrg(label="NODE", toler=tol)
            mapdl.allsel()

            # Z normal plane interface
            interface_z = i*Lz
            mapdl.nsel("S", "LOC", "Z", interface_z)
            mapdl.nummrg(label="NODE", toler=tol)
            mapdl.allsel()

        ##----- SYMMETRY GENERATION
        print("Generating symmetry")

        mapdl.allsel() 
        mapdl.nsel("S", "LOC", "Z", dist_z) # select nodes from areas normal to Z at coordenate z=dist_z
        mapdl.cm("Disp_Z", "NODE") # save these nodes as 'symm_Z'
        mapdl.allsel()

        mapdl.nsel("S", "LOC", "Y", dist_y)
        mapdl.cm("Disp_Y", "NODE")
        mapdl.allsel()

        mapdl.nsel("S", "LOC", "X", dist_x) 
        mapdl.cm("Disp_X", "NODE") 
        mapdl.allsel() 
        print("Symmetry generated")

        mapdl.allsel()
        mapdl.nsel("S", "LOC", "Z", 0) # select nodes from areas normal to Z at coordenate z=0
        mapdl.cm("CC_Z", "NODE") # save these nodes as 'CC_Z'
        A_z = dist_x*dist_y
        mapdl.allsel()
        mapdl.nsel("S", "LOC", "X", 0)
        mapdl.cm("CC_X", "NODE")
        A_x = dist_y*dist_z
        mapdl.allsel()
        mapdl.nsel("S", "LOC", "Y", 0)
        mapdl.cm("CC_Y", "NODE")
        A_y = dist_z*dist_x
        mapdl.allsel()

        ##----- MESH
        ele_size = Lx/10 #element size using cearactheristic length

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
        delta_L = Lz*0.01 #1% cell length
        if dir == "X": # For displacement in X, plane normal to Z (XZ)
            solution_G_modif(mapdl, "CC_Z", "Disp_Z", f"U{dir}", delta_L, 0)
        if dir == "Z": # For displacement in Z, plane normal to Y (YZ)
            solution_G_modif(mapdl, "CC_Y", "Disp_Y", f"U{dir}", delta_L, 0)
        if dir == "Y": # For displacement in Y, plane normal to X (XY)
            solution_G_modif(mapdl, "CC_X", "Disp_X", f"U{dir}", delta_L, 0)
        ######################################################################
        ##                 /POST1 o /POST26 |POST-PROCESS|                  ##
        ######################################################################
        mapdl.post1()

        ##----- RESULTS
        n_sets = mapdl.post_processing.nsets # number of substeps

        mapdl.set(lstep="LAST", sbstep="LAST") # analyse last results
        mapdl.allsel()

        # # Deformed structure
        # mapdl.run("/RGB,INDEX,100,100,100,0")   # white background
        # mapdl.run("/RGB,INDEX,0,0,0,15")        # black text
        # mapdl.run("/REPLOT")                    
        # mapdl.view(1, 0,1,0)
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
        # elem_nodes = np.array(mapdl.mesh.elem)[:,10:14] # elements nodes
        # elems = np.zeros((len(elem_ID),5)) # matrix = [ID N1 N2 N3 N4 ...] for the HDF5 file
        elems [:,0] = elem_ID.T
        elems [:,1:] = elem_nodes
        mapdl.finish() # Close post-processing stage
        
        ######################################################################
        ##                             CALCULUS                             ##
        ######################################################################
        if dir == "X": # G_XZ
            g, uat, ammag, fff = prop_mec_G_modif(mapdl, delta_L, dist_z, dist_x*dist_y, "Z", dir)
        if dir == "Z": # G_YZ
            g, uat, ammag, fff = prop_mec_G_modif(mapdl, delta_L, dist_y, dist_x*dist_z, "Y", dir)
        if dir == "Y": # G_XY
            g, uat, ammag, fff = prop_mec_G_modif(mapdl, delta_L, dist_x, dist_y*dist_z, "X", dir)
        
        G.extend([g/1e9]) # Pa
        #close and clear all
        mapdl.finish()
        mapdl.exit()
        shutil.rmtree(working_dir)

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
