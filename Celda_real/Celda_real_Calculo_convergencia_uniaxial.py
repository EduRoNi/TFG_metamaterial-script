'''
PyMAPDL script to determine the effective mechanical properties of metamaterial RVEs via tensile testing in the three principal directions (E and nu).V1
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
path = os.getcwd() # curretn working directory path 
if len(sys.argv) >= 3:
    file_name = sys.argv[1]
    reports_dir = sys.argv[2]
    temporal = sys.argv[3]

archivo = os.path.abspath(f"{file_name}.x_t") # find file

tolerancia = 5e-3 #tolerance in GPa
diff = 1
structure = 6 # initialise with 3x3x3

evo_EX = [] # EX evolution
evo_EY = [] # EY evolution
evo_EZ = [] # EZ evolution
num_celdas = [] # number of cells evolution
intento = 0 

try: 
    while diff > tolerancia:
        intento += 1
        print(f"Try {intento}")
       
        Direcciones = ["X", "Y", "Z"]
        for dir in Direcciones:           
            working_dir = os.path.join(temporal, f"CeldaReal.intent{intento}_{dir}") # Create a temporary directory to store and delete obsolete files
            mapdl = launch_mapdl(run_location=working_dir, override=True)
            
            mapdl.finish()
            mapdl.clear()

            ######################################################################
            ##                          INITIALISATION                          ##
            ######################################################################
            file_results = rf"{path}\{file_name}_results_{structure}x{structure}x{structure}"
            file_geometry = rf"{path}\{file_name}_geometry_{structure}x{structure}x{structure}"
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
            if intento == 1:
                rho_eff = rho(mapdl, Lx, Ly, Lz, 4430)
                print(f"Effective rho = {rho_eff} kg/m^3")

            # Delete and clean everything as we already have all relevant geometry data. Now open model previously meshed
            mapdl.finish()
            mapdl.clear()
            mapdl.cdread('db', 'ds.cdb')

            mapdl.slashsolu() # open solver
            mapdl.allsel()
            mapdl.ddele("ALL", "ALL") # delete all boundary conditions
            mapdl.finish()

            mapdl.prep7() # open pre-processor
            print(f"Elements: {mapdl.mesh.n_elem}\nNodes: {mapdl.mesh.n_node}")

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
            mapdl.cm("symm_Z", "NODE") # save these nodes as 'symm_Z'
            mapdl.allsel()

            mapdl.nsel("S", "LOC", "Y", dist_y)
            mapdl.cm("symm_Y", "NODE")
            mapdl.allsel()

            mapdl.nsel("S", "LOC", "X", dist_x) 
            mapdl.cm("symm_X", "NODE") 
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
            delta_L = Lz*0.001 # 1% cell length
            solution_1D_modif(mapdl, f"CC_{dir}", f"U{dir}", delta_L, 0)

            ######################################################################
            ##                 /POST1 o /POST26 |POST-PROCESS|                  ##
            ######################################################################
            mapdl.post1()

            ##----- RESULTS
            n_sets = mapdl.post_processing.nsets # number of substeps

            mapdl.set(lstep="LAST", sbstep="LAST") # analyse last results
            mapdl.allsel()

            # Deformed stucture
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
            # elem_nodes = np.array(mapdl.mesh.elem)[:,10:14] # elements nodes
            # elems = np.zeros((len(elem_ID),5)) # matrix = [ID N1 N2 N3 N4 ...] for the HDF5 file
            elems [:,0] = elem_ID.T
            elems [:,1:] = elem_nodes
            
            ######################################################################
            ##                             CALCULUS                             ##
            ######################################################################
            
            if dir == "X":
                L = dist_x
                E,s,e, nu_x, fff= prop_mec_modif(mapdl, delta_L, L, A_x, dir) # Young's Modulus X direction
                evo_EX.extend([E/1e9])
            if dir == "Y":
                L = dist_y
                E,s,e, nu_y, fff= prop_mec_modif(mapdl, delta_L, L, A_y, dir) # Young's Modulus Y direction
                evo_EY.extend([E/1e9])
            if dir == "Z":
                L = dist_z
                E,s,e, nu_z, fff= prop_mec_modif(mapdl, delta_L, L, A_z, dir) # Young's Modulus Z direction
                evo_EZ.extend([E/1e9])
            
            mapdl.finish() # Close post-processing stage

            # Close session and clean up temporary files
            mapdl.exit()
            shutil.rmtree(working_dir)
        
        
        # Check convergence tolerance
        if intento != 1:  
            diff = abs(E-E_viejo)/1e9 
            print(f"diff={diff}")
                
        if intento == 15:
            diff=tolerancia-1
        
        if diff > tolerancia:    
            structure = structure + 2 
            E_viejo = E
           
        num_celdas.extend([n])


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
    reports_dir_cell = os.path.join(reports_dir, file_name)
    if os.path.exists(reports_dir_cell) == False:
        os.makedirs(reports_dir_cell) # create directory to save final reports, in case it does not exist
    shutil.move(f"reporte_{file_name}.xlsx", os.path.join(reports_dir_cell, f"reporte_{file_name}.xlsx")) # move the file to correct directory


    print("File XLSX sucessfully saved.")
    print(f"Total run time={fin_programa-inicio_programa}")
    print(f"EX = {evo_EX} \nEY = {evo_EY} \nEZ = {evo_EZ}")