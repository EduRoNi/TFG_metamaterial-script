'''
PyMAPDL script to determine the elastoplastic behavior of metamaterial RVEs via tensile testing in the three principal directions.
Adapted to CEldillaV1 real.
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

# file_name = "CeldillaV1_real_II"
archivo = os.path.abspath(f"{file_name}.x_t") # find the file
# reports_dir = r"D:\Eduardo_Rodriguez\Python_archivos\Bucle\Reportes" # save all final reports for each cell
reports_dir_cell = os.path.join(reports_dir, file_name)
file_path = os.path.join(reports_dir_cell, f"reporte_{file_name}.xlsx")

# Read file to extract infromation about the cell
meta = pd.read_excel(file_path, sheet_name="Info")
structure = meta.loc[meta["Parámetro"] == "Final structure cells", "Valor"].values[0] #complete structure (later reduced to symmetries)
Ex = meta.loc[meta["Parámetro"] == "EX (GPa)", "Valor"].values[0]
Ey = meta.loc[meta["Parámetro"] == "EY (GPa)", "Valor"].values[0]
Ez = meta.loc[meta["Parámetro"] == "EZ (GPa)", "Valor"].values[0]
E_array = [Ex*1e9, Ey*1e9, Ez*1e9 ]  # E_array = [EX, EY, EZ] [Pa]
print(E_array)

direcciones = ["X", "Y", "Z"]
for j in range(0,len(direcciones)):
    dir = direcciones[j]
    print(f"ANALYSE IN {dir}")
    working_dir = os.path.join(temporal, f"CeldaReal.plast_{dir}") # Create a temporary directory to store and delete obsolete files
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
    mapdl.cm("Symm_Z", "NODE") # save these nodes as 'symm_Z'
    mapdl.allsel()

    mapdl.nsel("S", "LOC", "Y", dist_y)
    mapdl.cm("Symm_Y", "NODE")
    mapdl.allsel()

    mapdl.nsel("S", "LOC", "X", dist_x) 
    mapdl.cm("Symm_X", "NODE") 
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
    mapdl.emodif("ALL", "MAT", 5)

    ##----- SAVE GEOMETRY
    # mapdl.igesout(fname=file_geometry, ext="igs", att="1")
    # print(f"Geometry file saved as: {file_geometry}")

    mapdl.allsel()
    mapdl.finish() # Close pre-processing stage

    
    ######################################################################
    ##            /SOLU |SOLVER| y /POST1 |POSTPROCESS|                 ##
    ######################################################################
    if dir == "X":
        A = A_x # área m^2
        dist = dist_x # longitud característica en esa dirección
    if dir == "Y":
        A = A_y # área m^2
        dist = dist_y # longitud característica en esa dirección
    if dir == "Z":
        A = A_z # área m^2
        dist = dist_z # longitud característica en esa dirección

    rotura = False # Define solver stopping variable
    increment = dist*0.01 # 1% structure's characteristic length
    delta_L = 0 # imposed displacement
    mapdl.filname() # Changes the Jobname for the analysis to "work_dir"
    mapdl.slashsolu()
    mapdl.nlgeom("ON")
    mapdl.eqslv("SPARSE") 
    mapdl.nropt("FULL") # use FULL Newton-Rapshon method

    # Arrays for saving data
    sigma_global = [] # global sigma values
    epsilon_global = [] # global epsilon values
    FR = [] # reaction force
    disp_global = [] # global displacement
    load_step = 0 # load steps done
    old_substp = 0 # total substeps

    primera_vez = True

    # Symmetry boundary conditions
    mapdl.d("symm_X", "UX", 0)
    mapdl.d("symm_Y", "UY", 0)
    mapdl.d("symm_Z", "UZ", 0)

    while  rotura == False:
        try:
            mapdl.allsel()
            delta_L += increment # increment of imposed displacement
            mapdl.cmsel("S", f"CC_{dir}", "NODE")
            mapdl.d(node="ALL", lab=f"U{dir}", value=-delta_L)

            if primera_vez == True: # First iteration, initilise simulation
                print("new")
                mapdl.antype("STATIC", "NEW")
                mapdl.rescontrol("", "ALL", 1) # Save file for restart secuence
                primera_vez = False
            else: # Next iteration, restart from the previous load step
                print("restart")
                print(load_step, new_substp)
                mapdl.antype("", "REST", load_step, new_substp) # start from the last substep of the last step
                mapdl.cmsel("S", f"CC_{dir}", "NODE")
                mapdl.d(node="ALL", lab=f"U{dir}", value=-delta_L)

            mapdl.nsubst(10,20,5) # substeps definition
            mapdl.autots("ON")  # automatic time stepping ON 
            mapdl.outres("ALL", "ALL") #save all substeps

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

            mapdl.cmsel("S", f"CC_{dir}", "NODE")
            nodo_ref = mapdl.mesh.nnum[1] # reference node
            load_step += 1 
            new_substp = n_sets - old_substp
            old_substp = n_sets  
            for i in range(1, new_substp+1):
                mapdl.finish()
                mapdl.post1()
                mapdl.set(lstep= load_step, sbstep=str(i))
                mapdl.allsel()
                displacement,stress,strain,plastic=results(mapdl, node_ID, elem_ID)
                sigma_local = max(stress[:,7])
                plastic_local = max(plastic[:,1])
                ##########--------------------------################----------------------###########
                eppl = plastic[:,1] # equivalent plastic strain
                eppl_acc = plastic[:,2] # accumulated plastic strain
                plastic_values = eppl[eppl > plastic_lim]
                plastic_values_acc = eppl_acc[eppl_acc > plastic_lim] 
                print(len(plastic_values))
                print(len(plastic_values_acc))
                ##########--------------------------################----------------------###########
                pos = np.where(displacement[:,0] == nodo_ref)[0]
                delta = abs(displacement[int(pos[0]), j+1]) # actual displacement

                EZ, s, e, nuuuu, f= prop_mec_modif(mapdl, delta, dist, A, dir)
                print(f"For a total displacement of {delta} m, a local sigma = {sigma_local} and plastic strain = {plastic_local} are achieved. E{dir} = {EZ/1e9} GPa")
                sigma_global.extend([s])
                epsilon_global.extend([e])
                FR.extend([f])
                disp_global.extend([delta])
                if len(plastic_values_acc) > len(elem_ID) * 0.05:
                    print(f"The plastic strain limit has been exceeded (plastic limit = {plastic_lim}) in {len(plastic_values)} elements exceeding the 5% of elements damage with a global displacement of {delta} m (theoretically {delta_L} m).")
                    rotura = True
                    break
            if len(plastic_values) > len(elem_ID) * 0.05:
                break

            mapdl.finish()
            mapdl.slashsolu() # return to SOLVER for new iteration
        except:
            print("A problem has occurred")
            rotura = True
    
    print("New simulaitions")
    mapdl.finish()

    sigma_global = np.array(sigma_global) 
    epsilon_global = np.array(epsilon_global) 
    FR = np.array(FR)
    disp_global = np.array(disp_global)

    E = E_array[j] # [Pa]
    S, inte= plast(E, sigma_global, epsilon_global)
    print(f"Young's Modulus = {E/1e9} GPa, Yield stress S_{dir} = {S/1e6} MPa" )
    print(f"Max internal stress = {max(stress[:,7])/1e6} MPa")

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

    ax.plot(epsilon_global, sigma_global/1e6, color='blue', linewidth=2, label='Stress-Strain') 
    
    ax.scatter(epsilon_global[inte], S/1e6, marker='o', facecolors='none', edgecolors='black', 
                s=60, linewidths=1.2, zorder=5,
                label=f"Yield stress {np.round(S/1e6, 2)} MPa") 
    
    ax.scatter(epsilon_global[-1], sigma_global[-1]/1e6, marker='x', color='black', 
                s=60, linewidths=1.5, zorder=5,
                label=f"Failure {np.round(sigma_global[-1]/1e6, 2)} MPa") 
    
    ax.set_title(f'Stress-Strain in {dir}', fontsize=14, pad=10)
    ax.set_xlabel('Strain ($\epsilon$)', fontsize=13, labelpad=8)
    ax.set_ylabel('$\sigma$ [MPa]', fontsize=13, labelpad=8)
    
    ax.tick_params(top=True, right=True, direction='in')
    
    ax.legend(loc='best', edgecolor='black', framealpha=1)
    
    ax.grid(True, linestyle='-', linewidth=0.5, color='lightgray')
    
    plt.tight_layout()
    
    plt.savefig(f"Curva_TensionDeformacion_{dir}.png", dpi=300)
    plt.close()
    shutil.move(f"Curva_TensionDeformacion_{dir}.png", os.path.join(reports_dir_cell, f"Curva_TensionDeformacion_{dir}.png")) # move to correct directory


    # force-dsplacement
    fig, ax = plt.subplots(figsize=(7, 6))

    ax.plot(disp_global, FR/1e3, color='blue', linewidth=2, label='Force-Displacement') 

    ax.scatter(disp_global[-1], FR[-1]/1e3, marker='x', color='black', 
                s=60, linewidths=1.5, zorder=5,
                label=f"Failure {np.round(FR[-1]/1e3, 2)} kN") 
    
    ax.set_title(f'Force-Displacement in {dir}', fontsize=14, pad=10)
    ax.set_xlabel('Displacement [mm]', fontsize=13, labelpad=8)
    ax.set_ylabel('Force [kN]', fontsize=13, labelpad=8)
    
    ax.tick_params(top=True, right=True, direction='in')
    
    ax.legend(loc='best', edgecolor='black', framealpha=1)
    
    ax.grid(True, linestyle='-', linewidth=0.5, color='lightgray')
    
    plt.tight_layout()
    
    plt.savefig(f"Curva_FuerzaDesplaz_{dir}.png", dpi=300)
    plt.close()
    shutil.move(f"Curva_FuerzaDesplaz_{dir}.png", os.path.join(reports_dir_cell, f"Curva_FuerzaDesplaz_{dir}.png")) # move to correct directory

    # File control
    # Meta
    meta = pd.read_excel(file_path, sheet_name="Info")
    meta_new = pd.concat([meta, pd.DataFrame({
        "Parámetro": [f"Límite elast_{dir} (MPa)", f"Tensión de rotura_{dir} (MPa)", f"Tensión interna máx (MPa)"],
        "Valor": [S/1e6, sigma_global[-1]/1e6, sigma_local/1e6]
    })], ignore_index=True)

    # Overwrite to include new information
    with pd.ExcelWriter(file_path, engine='openpyxl', mode='a', if_sheet_exists='replace') as writer:
        meta_new.to_excel(writer, sheet_name="Info", index=False)

        # Data
        datos_plast = {
            f"Sigma_{dir} (Pa)": sigma_global,
            f"Epsilon_{dir}": epsilon_global, 
            f"Displacement_{dir} (mm)": disp_global, 
            f"Force_{dir} (N)": FR
        }
        df_plast = pd.DataFrame(datos_plast)
        df_plast.to_excel(writer, sheet_name=f"Plasticidad_{dir}", index=False)

    print("File XLSX sucessfully saved.")

    print(f"FINISHED ANALYSIS IN {dir}")
    mapdl.finish()
    mapdl.exit()
    shutil.rmtree(working_dir)

# except Exception as e: # if any error occurs
#     print(e)
#     mapdl.finish()
#     mapdl.exit()
#     shutil.rmtree(working_dir)

# finally:
#     fin_programa = time.time()
#     print(f"Total time={fin_programa-inicio_programa}")