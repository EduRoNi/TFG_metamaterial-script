"""
This file collects all the functions used across the various Python files, as well as all the versions created for specific use cases. 
"""

import h5py 
import numpy as np


# Extract results. Used in every file.
def results(mapdl, node_ID, elem_ID):
    '''
    Function that extracts simulation data for a specific substep.
    
    INPUTS
    node_ID: vector containing all node IDs using mapdl.mesh.nnum
    elem_ID: vector containing all element IDs using mapdl.mesh.enum

    OUTPUTS
    displacement: displacement of each node
    stress: stresses at the element centroids
    strain: equivalent elastic strains at the element centroids
    plastic: equivalent plastic strains at the element centroids
    # teqv: Total Equivalent Strain, summing elastic+plastic+thermal+creep+swelling
    '''
    # Nodal displacements [m]
    u = mapdl.post_processing.nodal_displacement("ALL")
    u_norm = mapdl.post_processing.nodal_displacement("NORM")
    displacement = np.zeros((len(node_ID), 5)) # matrix = [ID ux uy uz norm]
    displacement [:,0] = node_ID.T
    displacement [:,1:4] = u
    displacement [:,4] = u_norm.T

    # Elements stresses [Pa]
    # 'AVG' store averaged element *centroid value* of the specified item component
    sx  = mapdl.post_processing.element_stress("X", "AVG")
    sy  = mapdl.post_processing.element_stress("Y", "AVG")
    sz  = mapdl.post_processing.element_stress("Z", "AVG")
    sxy = mapdl.post_processing.element_stress("XY", "AVG")
    syz = mapdl.post_processing.element_stress("YZ", "AVG")
    sxz = mapdl.post_processing.element_stress("XZ", "AVG")
    seqv = mapdl.post_processing.element_stress("EQV", "AVG")  # Von Mises
    stress = np.zeros((len(elem_ID), 8)) # matrix = [ID sx sy sz sxy syz sxz seqv]
    stress [:,0] = elem_ID.T
    stress [:,1] = sx.T
    stress [:,2] = sy.T
    stress [:,3] = sz.T
    stress [:,4] = sxy.T
    stress [:,5] = syz.T
    stress [:,6] = sxz.T
    stress [:,7] = seqv.T

    # Elements elastic strains
    # 'AVG' store averaged element *centroid value* of the specified item component
    ex  = mapdl.post_processing.element_values("EPEL", "X", "AVG")
    ey  = mapdl.post_processing.element_values("EPEL", "Y", "AVG")
    ez  = mapdl.post_processing.element_values("EPEL", "Z", "AVG")
    exy = mapdl.post_processing.element_values("EPEL", "XY", "AVG")
    eyz = mapdl.post_processing.element_values("EPEL", "YZ", "AVG")
    exz = mapdl.post_processing.element_values("EPEL", "XZ", "AVG")
    eeqv = mapdl.post_processing.element_values("EPEL", "EQV") # Elastic equivalent strain
    strain = np.zeros((len(elem_ID), 8)) # matrix = [ID ex ey ez exy eyz exz eeqv]
    strain [:,0] = elem_ID.T
    strain [:,1] = ex.T
    strain [:,2] = ey.T
    strain [:,3] = ez.T
    strain [:,4] = exy.T
    strain [:,5] = eyz.T
    strain [:,6] = exz.T
    strain [:,7] = eeqv.T

    # Elements plastic strains
    plc = mapdl.post_processing.element_values("EPPL", "EQV", "AVG") # equivalent plastic strain
    plc_acc = mapdl.post_processing.element_values("NL", "EPEQ", "AVG") # accumulated plastic strain
    plastic = np.zeros((len(elem_ID), 3))
    plastic [:,0] = elem_ID.T
    plastic [:,1] = plc.T
    plastic [:,2] = plc_acc.T

    # # Equivalent total strains
    # teqv = mapdl.post_processing.element_values("EPTT", "EQV") # Total equivalent strain

    return displacement, stress, strain, plastic

# Extend de 'results_file' for AI training. Used in every file.
def escritura_h5(name, displacement, stress, strain, plastic, sbtep):
    '''
    Function that opens the subsequently created h5 document to add the data from each substep

    INPUTS
    name: file name
    displacement: displacement of each node
    stress: element stresses
    strain: element elastic strains
    plastic: element equivalent plastic strains
    step: substep number being analyzed

    OUTPUTS
    '''
    name_sb = f"Step {sbtep}"
    with h5py.File(name, "a", track_order=True) as f:
        # Substep
        stp = f.create_group(name_sb)
        stp.create_dataset("Displacements", data=displacement)
        stp.create_dataset("Stresses", data=stress)
        stp.create_dataset("Elastic Strain", data=strain)
        stp.create_dataset("Plastic Strain", data=plastic)

# Determine elastic properties, E and nu. 
def prop_mec(mapdl, delta_L, L, A, dirección):
    '''
    Function that calculates the mechanical properties of the structure after a tensile test

    
    INPUTS
    delta_L: prescribed displacement [m]
    L: length of the structure in the direction of the prescribed displacement [m]
    A: value of the area affected by the prescribed displacement. It must be the area considering a solid block, including the voids. [m^2]
    direction: principal direction in which the property is calculated. Valid options: "X", "Y" or "Z" (string)

    OUTPUTS
    E: Young s modulus [Pa]
    epsilon: global strain (macro)
    sigma: global stress [Pa] (macro)
    '''
    # Young's Modulus (E)
    mapdl.cmsel("S", f"CC_{dirección}", "AREA")
    mapdl.nsla("S", 1)
    mapdl.post1()
    mapdl.fsum() #Sums the nodal force and moment contributions of elements. Sum all nodal forces in global Cartesian coordinate system (default). Sum all nodal forces for all selected nodes (default), excluding contact elements.
    F = mapdl.get_value(entity="FSUM", item1="ITEM", it1num=f"F{dirección}") # reaction force
    mapdl.finish()
    epsilon = delta_L/L
    sigma = F / A
    E = abs(sigma / epsilon)
    print(f"E={E/1e9} GPa")
    # print(f"epsilon={epsilon}, sigma={sigma/1e6} MPa, F={F} N, A={A} m^2")

    # Poisson coeficient (nu)
    if dirección == "Z":
        dir2 = ["X", "Y"]
    elif dirección == "X":
        dir2 = ["Y" , "Z"]
    elif dirección == "Y":
        dir2 = ["X", "Z"]
    nu = []
    for dir in dir2:
        mapdl.allsel() 
        mapdl.cmsel("S", f"CC_{dir}", "AREA")
        mapdl.nsla("S", 1)
        u_trans = 0
        for node in mapdl.mesh.nnum:
            u_trans += mapdl.get_value("NODE", node, "U", dir)
        u_trans = - u_trans/len(mapdl.mesh.nnum)
        epsilon_trans = u_trans/L
        nu.extend([-epsilon_trans/epsilon])
    print(f"nu_{dirección}{dir}={nu}")    
    # print(f"u_trans={u_trans} m, epsilon_trans={epsilon_trans}")
    
    return E, abs(sigma), abs(epsilon), np.array(nu), abs(F)
def prop_mec_modif(mapdl, delta_L, L, A, dirección): # Adapted for CeldillaV1_real
    '''
    Function that calculates the mechanical properties of the structure after a tensile test.
    Modified and adapted for CeldillaV1_real.

    INPUTS
    delta_L: prescribed displacement [m]
    L: length of the structure in the direction of the prescribed displacement [m]
    A: value of the area affected by the prescribed displacement. It must be the area considering a solid block, including the voids. [m^2]
    direction: principal direction in which the property is calculated. Valid options: "X", "Y" or "Z" (string)

    OUTPUTS
    E: Young s modulus [Pa]
    epsilon: global strain (macro)
    sigma: global stress [Pa] (macro)
    '''
    # Young's Modulus (E)
    mapdl.cmsel("S", f"CC_{dirección}", "NODE")
    mapdl.post1()
    mapdl.fsum() #Sums the nodal force and moment contributions of elements. Sum all nodal forces in global Cartesian coordinate system (default). Sum all nodal forces for all selected nodes (default), excluding contact elements.
    F = mapdl.get_value(entity="FSUM", item1="ITEM", it1num=f"F{dirección}")
    mapdl.finish()
    epsilon = delta_L/L
    sigma = F / A
    E = abs(sigma / epsilon)
    print(f"E={E/1e9} GPa")
    # print(f"epsilon={epsilon}, sigma={sigma/1e6} MPa, F={F} N, A={A} m^2")

    # Poisson coeficient (nu)
    if dirección == "Z":
        dir2 = ["X", "Y"]
    elif dirección == "X":
        dir2 = ["Y" , "Z"]
    elif dirección == "Y":
        dir2 = ["X", "Z"]
    nu = []
    for dir in dir2:
        mapdl.allsel() 
        mapdl.cmsel("S", f"CC_{dir}", "NODE")
        u_trans = 0
        for node in mapdl.mesh.nnum:
            u_trans += mapdl.get_value("NODE", node, "U", dir)
        u_trans = - u_trans/len(mapdl.mesh.nnum)
        epsilon_trans = u_trans/L
        nu.extend([-epsilon_trans/epsilon])
    print(f"nu_{dirección}{dir}={nu}")    
    # print(f"u_trans={u_trans} m, epsilon_trans={epsilon_trans}")
    
    return E, abs(sigma), abs(epsilon), np.array(nu), abs(F)

# Determine elastic properties, G. Last version in Metamaterial_G.
def prop_mec_G(mapdl, delta, L, A, dir1, dir2):
    '''
    Calculates the shear modulus G [Pa]

    INPUTS
    delta: prescribed displacement [m]
    L: dimension transversal to the prescribed displacement [m]
    A: area of the surface where the displacement is prescribed [m^2]
    dir1: direction perpendicular to the face where the displacement is prescribed (L)
    dir2: direction in which the displacement is prescribed (delta)

    OUTPUTS
    G: shear modulus [Pa]
    '''
    # Surface nodes
    mapdl.cmsel("S", f"Disp_{dir1}", "AREA")
    mapdl.nsla("S", 1)
    F = 0
    for node in mapdl.mesh.nnum:
        F += mapdl.get_value("NODE", node, "RF", f"F{dir2}") # reaction force
    gamma = delta/L
    tau = F / A
    G = abs(tau / gamma)
    print(f"G_{dir1}{dir2}= {G/1e9} GPa")

    return G, tau, gamma, abs(F)
def prop_mec_G_modif(mapdl, delta, L, A, dir1, dir2): # Adapted for CeldillaV1_real
    '''
    Calculates the shear modulus G [Pa] 
    Modified and adapted for CeldillaV1_real.  
    
    INPUTS
    delta: prescribed displacement [m]
    L: dimension transversal to the prescribed displacement [m]
    A: area of the surface where the displacement is prescribed [m^2]
    dir1: direction perpendicular to the face where the displacement is prescribed (L)
    dir2: direction in which the displacement is prescribed (delta)

    OUTPUTS
    G: shear modulus [Pa]
    '''
    # Surface nodes
    mapdl.cmsel("S", f"Disp_{dir1}", "NODE")
    F = 0
    for node in mapdl.mesh.nnum:
        F += mapdl.get_value("NODE", node, "RF", f"F{dir2}") # reaction force
    gamma = delta/L
    tau = F / A
    G = abs(tau / gamma)
    print(f"G_{dir1}{dir2}= {G/1e9} GPa")

    return G, tau, gamma, abs(F)

# Solver para cunado resulevo en las tres direcciones usando el mismo launch (incluye .ddele()). Usado en Metamaterial_calculo_en_las_3_direcciones.
def solution_1D(mapdl, superficie, direccion, delta_L, sol):
    '''
    Function that calls the APDL /SOLU processor and solves load cases with prescribed displacement on surfaces in **one direction**.

    INPUTS
    superficie: surface where the displacement is prescribed. Valid options: "CC_X", "CC_Y", "CC_Z" (the ones saved in /PREP)
    direccion: direction in which the displacement is prescribed. Valid options: "UX", "UY", "UZ", "ROTX", "ROTY", "ROTZ" (the ones allowed by the MAPDL .d() module)
    delta_L: value of the prescribed displacement
    sol: linear or non-linear solution type with substeps. 0-Linear. 1-Non-linear (substep). 2-Fixed  
    '''
    mapdl.slashsolu() # /SOLU processor
    mapdl.allsel()
    mapdl.ddele("ALL", "ALL") # deletes boundary conditions

    mapdl.cmsel("S", superficie, "AREA")
    mapdl.nsla("S", 1) 
    mapdl.cm("CC", "NODE")
    mapdl.d(node="ALL", lab=direccion, value=-delta_L)
    print(f"displacement: {-delta_L}")

    ##----- SOLVER CONTROLS
    # SOLVER CONTROLS FOR LINEAL SOLUTIONS
    if sol == 0:
        print("Lineal Solution")
        mapdl.antype(antype="STATIC", status="NEW")
        mapdl.nlgeom(key="OFF") #large deflection OFF
        # mapdl.nropt()
        mapdl.eqslv("SPARSE") 
        mapdl.nsubst(1,1,1) #substeps definition 
        mapdl.autots("OFF")  # automatic time stepping OFF
        ##----- SOLVER
        print("Solver initialised")
        mapdl.allsel()
        mapdl.solve()
        print("Solver finished")

    # SOLVER CONTROLS FOR NONLINEAR SOLUTION
    if sol == 1:
    #Load simluation
        print("Non Lineal Solution")
        mapdl.antype("STATIC", "NEW")
        mapdl.nlgeom("ON") #large deflection ON
        mapdl.eqslv("SPARSE") 
        mapdl.nropt("FULL") #use FULL Newton-Rapshon method
        mapdl.autots("OFF") #automatic time stepping OFF
        mapdl.nsubst(50, 100, 25) #substep definition
        mapdl.cnvtol("U", toler=0.05, norm=2)  # Displacement tolerance 5%, the one recommended by the guide for displacements
        mapdl.outres("ALL", "ALL") #save all substeps
        ##----- SOLVER
        print("Solver initialised")
        mapdl.allsel()
        mapdl.solve()
        print("Solver finished")

    mapdl.finish()
def solution_1D_modif(mapdl, superficie, direccion, delta_L, sol): # Adapted for CeldillaV1_real
    '''
    Function that calls the APDL /SOLU processor and solves load cases with prescribed displacement on surfaces in **one direction**.
    Modified and adapted for CeldillaV1_real.

    INPUTS
    superficie: surface where the displacement is prescribed. Valid options: "CC_X", "CC_Y", "CC_Z" (the ones saved in /PREP)
    direccion: direction in which the displacement is prescribed. Valid options: "UX", "UY", "UZ", "ROTX", "ROTY", "ROTZ" (the ones allowed by the MAPDL .d() module)
    delta_L: value of the prescribed displacement
    sol: linear or non-linear solution type with substeps. 0-Linear. 1-Non-linear (substep). 2-Fixed  
    '''
    mapdl.slashsolu() # /SOLU processor
    mapdl.allsel()
    mapdl.ddele("ALL", "ALL") # Delete boundary conditions

    # Symmetry boundary conditions
    mapdl.d("symm_X", "UX", 0)
    mapdl.d("symm_Y", "UY", 0)
    mapdl.d("symm_Z", "UZ", 0)

    mapdl.cmsel("S", superficie, "NODE")
    mapdl.cm("CC", "NODE")
    mapdl.d(node="ALL", lab=direccion, value=-delta_L)
    print(f"displacement: {-delta_L}")

    ##----- SOLVER CONTROLS
    # SOLVER CONTROLS FOR LINEAL SOLUTIONS
    if sol == 0:
        print("Lineal Solution")
        mapdl.antype(antype="STATIC", status="NEW")
        mapdl.nlgeom(key="OFF") #large deflection OFF
        # mapdl.nropt()
        mapdl.eqslv("PCG") # method
        mapdl.nsubst(1,1,1) #substeps definition 
        mapdl.autots("OFF")  # automatic time stepping OFF
        ##----- SOLVER
        print("Solver initialised")
        mapdl.allsel()
        mapdl.solve()
        print("Solver finished")

    # SOLVER CONTROLS FOR NONLINEAR SOLUTION
    if sol == 1:
    #Load simluation
        print("Non Lineal Solution")
        mapdl.antype("STATIC", "NEW")
        mapdl.nlgeom("ON") #large deflection ON
        mapdl.eqslv("SPARSE") 
        mapdl.nropt("FULL") #use FULL Newton-Rapshon method
        mapdl.autots("OFF") #automatic time stepping OFF
        mapdl.nsubst(50, 100, 25) #substep definition
        mapdl.cnvtol("U", toler=0.05, norm=2)  # Displacement tolerance 5%, the one recommended by the guide for displacements
        mapdl.outres("ALL", "ALL") #save all substeps
        ##----- SOLVER
        print("Solver initialised")
        mapdl.allsel()
        mapdl.solve()
        print("Solver finished")

    mapdl.finish()

# Solver para resolver en una sola dirección (no usa .ddele()). NO se usa
def solution(mapdl, superficie, direccion, delta_L, sol):
    '''
    Function that calls the APDL /SOLU processor and solves load cases with prescribed displacement on surfaces in **one direction** with no boundary conditions elimination.

    INPUTS
    superficie: surface where the displacement is prescribed. Valid options: "CC_X", "CC_Y", "CC_Z" (the ones saved in /PREP)
    direccion: direction in which the displacement is prescribed. Valid options: "UX", "UY", "UZ", "ROTX", "ROTY", "ROTZ" (the ones allowed by the MAPDL .d() module)
    delta_L: value of the prescribed displacement
    sol: linear or non-linear solution type with substeps. 0-Linear. 1-Non-linear (substep). 2-Fixed  
    '''
    mapdl.slashsolu() 
    mapdl.cmsel("S", superficie, "AREA")
    mapdl.nsla("S", 1) 
    mapdl.cm("CC", "NODE") 
    mapdl.d(node="ALL", lab=direccion, value=-delta_L)
    print(f"displacement: {-delta_L}")

    ##----- SOLVER CONTROLS
    # SOLVER CONTROLS FOR LINEAL SOLUTIONS
    if sol == 0:
        print("Lineal Solution")
        mapdl.antype(antype="STATIC", status="NEW")
        mapdl.nlgeom(key="OFF") #large deflection OFF
        # mapdl.nropt()
        mapdl.eqslv("SPARSE") 
        mapdl.nsubst(1,1,1) #substeps definition 
        mapdl.autots("OFF")  # automatic time stepping OFF
        ##----- SOLVER
        print("Solver initialised")
        mapdl.allsel()
        mapdl.solve()
        print("Solver finished")

    # SOLVER CONTROLS FOR NONLINEAR SOLUTION
    if sol == 1:
    #Load simluation
        print("Non Lineal Solution")
        mapdl.antype("STATIC", "NEW")
        mapdl.nlgeom("ON") #large deflection ON
        mapdl.eqslv("SPARSE") 
        mapdl.nropt("FULL") #use FULL Newton-Rapshon method
        mapdl.autots("OFF") #automatic time stepping OFF
        mapdl.nsubst(50, 100, 25) #substep definition
        mapdl.cnvtol("U", toler=0.05, norm=2)  # Displacement tolerance 5%, the one recommended by the guide for displacements
        mapdl.outres("ALL", "ALL") #save all substeps
        ##----- SOLVER
        print("Solver initialised")
        mapdl.allsel()
        mapdl.solve()
        print("Solver finished")

    mapdl.finish()

# Specific solver for simple shear conditions. Used in Metamaterial_G.
def solution_G(mapdl, sup_fixed, sup_disp, direccion, delta_L, sol):
    '''
    Function that calls the APDL /SOLU processor and solves load cases with prescribed displacement on surfaces for calculating G.
    Generates a pure shear condition for calculating G.

    INPUTS
    sup_fixed: surface where the FIXED condition is applied. Valid options: "CC_X", "CC_Y", "CC_Z" (the ones saved in /PREP)
    sup_disp: surface where the prescribed displacements are applied. Valid options: "Disp_Z", "Disp_X", "Disp_Y" (the ones saved in /PREP)
    direccion: direction in which the displacement is prescribed. Valid options: "UX", "UY", "UZ", "ROTX", "ROTY", "ROTZ" (the ones allowed by the MAPDL .d() module)
    delta_L: value of the prescribed displacement
    sol: linear or non-linear solution type with substeps. 0-Linear. 1-Non-linear (substep)
    '''
    mapdl.slashsolu() # /SOLU processor
    mapdl.allsel()
    mapdl.ddele("ALL", "ALL") # delete boundary conditions
    # FIXED
    mapdl.cmsel("S", sup_fixed, "AREA")
    mapdl.nsla("S", 1)
    mapdl.d(node="ALL", lab="ALL", value=0)

    # CONSTRAINTS LATERAL FACES
    if direccion == "UX":
        mapdl.cmsel("S", "Disp_X", "AREA")
        mapdl.cmsel("A", "CC_X", "AREA")
        mapdl.cmsel("A", "Disp_Y", "AREA")
        mapdl.cmsel("A", "CC_Y", "AREA")
        mapdl.nsla("S", 1)
        mapdl.d(node="ALL", lab="UY", lab2="UZ", value=0)
    if direccion == "UY":
        mapdl.cmsel("S", "Disp_Y", "AREA")
        mapdl.cmsel("A", "CC_Y", "AREA")
        mapdl.cmsel("A", "Disp_Z", "AREA")
        mapdl.cmsel("A", "CC_Z", "AREA")
        mapdl.nsla("S", 1)
        mapdl.d(node="ALL", lab2="UX", lab3="UZ", value=0)
    if direccion == "UZ":
        mapdl.cmsel("S", "Disp_Z", "AREA")
        mapdl.cmsel("A", "CC_Z", "AREA")
        mapdl.cmsel("A", "Disp_X", "AREA")
        mapdl.cmsel("A", "CC_X", "AREA")
        mapdl.nsla("S", 1)
        mapdl.d(node="ALL", lab2="UX", lab3="UY", value=0)

    # DISPLACEMENT
    mapdl.allsel()
    mapdl.cmsel("S", sup_disp, "AREA")
    mapdl.nsla("S", 1)
    if direccion == "UX":
        mapdl.d(node="ALL", lab=direccion, value=delta_L)
        mapdl.d(node="ALL", lab="UY", lab2="UZ", value=0)
    if direccion == "UY":
        mapdl.d(node="ALL", lab=direccion, value=delta_L)
        mapdl.d(node="ALL", lab2="UX", lab3="UZ", value=0)
    if direccion == "UZ":
        mapdl.d(node="ALL", lab=direccion, value=delta_L)
        mapdl.d(node="ALL", lab2="UX", lab3="UY", value=0)
    print(f"displacement: {delta_L}")

    ##----- SOLVER CONTROLS
    # SOLVER CONTROLS FOR LINEAL SOLUTIONS
    if sol == 0:
        print("Lineal Solution")
        mapdl.antype(antype="STATIC", status="NEW")
        mapdl.nlgeom(key="OFF") #large deflection OFF
        # mapdl.nropt()
        mapdl.eqslv("SPARSE") 
        mapdl.nsubst(1,1,1) #substeps definition 
        mapdl.autots("OFF")  # automatic time stepping OFF
        ##----- SOLVER
        print("Solver initialised")
        mapdl.allsel()
        mapdl.solve()
        print("Solver finished")

    # SOLVER CONTROLS FOR NONLINEAR SOLUTION
    if sol == 1:
    #Load simluation
        print("Non Lineal Solution")
        mapdl.antype("STATIC", "NEW")
        mapdl.nlgeom("ON") #large deflection ON
        mapdl.eqslv("SPARSE") 
        mapdl.nropt("FULL") #use FULL Newton-Rapshon method
        mapdl.autots("OFF") #automatic time stepping OFF
        mapdl.nsubst(50, 100, 25) #substep definition
        mapdl.cnvtol("U", toler=0.05, norm=2)  # Displacement tolerance 5%, the one recommended by the guide for displacements
        mapdl.outres("ALL", "ALL") #save all substeps
        ##----- SOLVER
        print("Solver initialised")
        mapdl.allsel()
        mapdl.solve()
        print("Solver finished")

    mapdl.finish()
def solution_G_modif(mapdl, sup_fixed, sup_disp, direccion, delta_L, sol): # Adapted for CeldillaV1_real
    '''
    Function that calls the APDL /SOLU processor and solves load cases with prescribed displacement on surfaces for calculating G.
    Generates a pure shear condition for calculating G.
    Modified and adapted for CeldillaV1_real.

    INPUTS
    sup_fixed: surface where the FIXED condition is applied. Valid options: "CC_X", "CC_Y", "CC_Z" (the ones saved in /PREP)
    sup_disp: surface where the prescribed displacements are applied. Valid options: "Disp_Z", "Disp_X", "Disp_Y" (the ones saved in /PREP)
    direccion: direction in which the displacement is prescribed. Valid options: "UX", "UY", "UZ", "ROTX", "ROTY", "ROTZ" (the ones allowed by the MAPDL .d() module)
    delta_L: value of the prescribed displacement
    sol: linear or non-linear solution type with substeps. 0-Linear. 1-Non-linear (substep)
    '''
    mapdl.slashsolu() 
    mapdl.allsel()
    mapdl.ddele("ALL")

    # FIXED
    mapdl.cmsel("S", sup_fixed, "NODE")
    mapdl.d(node="ALL", lab="ALL", value=0)
    # DISPLACEMENT
    mapdl.allsel()
    mapdl.cmsel("S", sup_disp, "NODE")
    if direccion == "UX":
        mapdl.d(node="ALL", lab=direccion, value=delta_L)
        mapdl.d(node="ALL", lab="UY", lab2="UZ", value=0)
    if direccion == "UY":
        mapdl.d(node="ALL", lab=direccion, value=delta_L)
        mapdl.d(node="ALL", lab2="UX", lab3="UZ", value=0)
    if direccion == "UZ":
        mapdl.d(node="ALL", lab=direccion, value=delta_L)
        mapdl.d(node="ALL", lab2="UX", lab3="UY", value=0)
    print(f"displacement: {+delta_L}")

    ##----- SOLVER CONTROLS
    # SOLVER CONTROLS FOR LINEAL SOLUTIONS
    if sol == 0:
        print("Lineal Solution")
        mapdl.antype(antype="STATIC", status="NEW")
        mapdl.nlgeom(key="OFF") #large deflection OFF
        # mapdl.nropt()
        mapdl.eqslv("SPARSE") 
        mapdl.nsubst(1,1,1) #substeps definition 
        mapdl.autots("OFF")  # automatic time stepping OFF
        ##----- SOLVER
        print("Solver initialised")
        mapdl.allsel()
        mapdl.solve()
        print("Solver finished")

    # SOLVER CONTROLS FOR NONLINEAR SOLUTION
    if sol == 1:
    #Load simluation
        print("Non Lineal Solution")
        mapdl.antype("STATIC", "NEW")
        mapdl.nlgeom("ON") #large deflection ON
        mapdl.eqslv("SPARSE") 
        mapdl.nropt("FULL") #use FULL Newton-Rapshon method
        mapdl.autots("OFF") #automatic time stepping OFF
        mapdl.nsubst(50, 100, 25) #substep definition
        mapdl.cnvtol("U", toler=0.05, norm=2)  # Displacement tolerance 5%, the one recommended by the guide for displacements
        mapdl.outres("ALL", "ALL") #save all substeps
        ##----- SOLVER
        print("Solver initialised")
        mapdl.allsel()
        mapdl.solve()
        print("Solver finished")

    mapdl.finish()

# Determine Sy for plasticity behavior. Used in Metamaterial_prueba_hasta_rotura_platicidad and Metamaterial_G_plasticidad_hasta_rotura
def plast(E, sigma, epsilon):
    '''
    Function to extract the plastic strain-stress pairs to feed the Multilinear Isotropic Hardening model.

    INPUTS
    E: Young s modulus of the structure (effective E, not the real material E) [Pa]
    sigma: array containing the global stress values under tension
    epsilon: total equivalent strain (including elastic and plastic components) from the tensile test

    OUTPUTS
    S_y: effective yield strength of the structure [Pa]
    data: plastic strain-stress pairs
    '''
    # Line with a 0.2% plastic strain offset
    eq_offset = E * (epsilon - 0.002)
    dif = sigma - eq_offset
    intersecciones = np.where(np.diff(np.sign(dif)))[0] # np.sign() devuelve -1, 1, 0 dependiendo si la resta da valor negativo, positivo o 0
                                                        # np.diff() clacula la diferencia entre un elemento y el anterior, de manera que ientras la curva 
                                                        #           y la recta no se crucen, los signos serán iguales (ej. 1, 1, 1), por lo que la resta será 0.
                                                        #           En el momento en que se cruzan, el signo cambia (de 1 a -1 o viceversa). La resta será 2 o -2.
                                                        # np.where() esta función nos da los índices (las posiciones en el array) donde el valor es distinto 
                                                        #           de cero. Es decir, nos dice en qué posición exacta de tus datos ocurrió el cruce.
    S_y = sigma[intersecciones[0]]
    
    # Data extraction for the multilinear model
    # The transition to plasticity occurs at the index "intersecciones[0]"
    def_plast = epsilon - sigma/E
    def_plast = def_plast[intersecciones[0]:] # stores the values from the intersection onwards
    data = np.zeros([len(def_plast),2]) # create the matrix to input the data pairs
    data [:,0] = def_plast.T 
    data [:,1] = sigma[intersecciones[0]:].T
    return S_y

# Determine effective rho of the cells. Used only in Matamaterial_calculo_en_las_3_direcciones.py
def rho(mapdl, Lx, Ly, Lz, rho):
    """
    Return the value of effective density of the cell

    INPUTS
    Lx, Ly, Lz: box dimensions of the cell [m]
    rho: density of the base material [kg/m^3]

    OUTPUTS
    rho_eff: effective density [kg/m^3]
    """
    mapdl.vsum()
    vol_real = mapdl.get_value("VOLU", 1, "VOLU")
    vol_eff = Lx*Ly*Lz
    rho_eff = (vol_real/vol_eff)*rho

    return rho_eff

# Determine the value of "n" from Gibson-Ashby model. Used only in Matamaterial_calculo_en_las_3_direcciones.py
def gibson_ashby(rho_base, E_base, rho_eff, E_eff):
    """
    Function that determine the value of "n" from de Gibson-Ashby model for lattice structures

    INPUTS
    rho_base: material density [kg/m^3]
    E_base: material Young s Modulus [GPa]
    rho_eff: effective density [kg/m^3]
    E_eff: cell effective Young s Modulus [GPa]

    OUTPUTS
    n: exponent from Gibson-Ashby model [-]
    """
    n = np.log(E_eff/E_base)/np.log(rho_eff/rho_base)
    return np.round(n, 3)

# Determine Jacobian Ratio and Aspect Ratio. Used only in elastic convergence codes.
def jacobian_ratio(mapdl):
    """
    Calculates the Jacobian Ratio of the model. It must be ran into PREP7 module, before any simulaition starts.
    
    INPUTS

    OUTPUTS
    jratio: Jacobian ratio
    """

    jratio= 0
    mapdl.allsel()
    elem_id = mapdl.mesh.enum
    for eid in elem_id:
        jac = mapdl.get_value("ELEM", eid, "SHPAR", "JACR")
        jratio += jac
    jratio = jratio/len(mapdl.mesh.enum)
    print(f"Jacobian Ratio = {jratio}")

    return jratio
def aspect_ratio(mapdl):
    """
    Calculates the Aspect Ratio of the model. It must be ran into PREP7 module, before any simulaition starts.
    
    INPUTS

    OUTPUTS
    aratio: Aspect ratio
    """
    aratio= 0
    mapdl.allsel()
    elem_id = mapdl.mesh.enum
    for eid in elem_id:
        aac = mapdl.get_value("ELEM", eid, "SHPAR", "ASPE")
        aratio += aac
    aratio = aratio/len(mapdl.mesh.enum)
    print(f"Aspect Ratio = {aratio}")
    
    return aratio