#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 10 15:45:21 2021

@author: alexhoffnagle
"""
import argparse
import numpy as np
from pyrosetta.rosetta.numeric import xyzVector_double_t

#Uses input parameters to define 3 coordinates that will be used to position the protein chain
def generate_coords_from_parameters(d1,d2,Angle1,Angle2,Angle3,Angle4,Angle5,resn):
    A1=np.array([0,0,0]) #Position atom 1 (His Nep or Asp/Glu OD2) at origin
    RAngle3x=np.array([[1,0,0],[0,np.cos(Angle3),-np.sin(Angle3)],[0,np.sin(Angle3),np.cos(Angle3)]]) #rotation matrix for rotation about x axis by Theta3 rad
    RAngle1z=np.array([[np.cos(Angle1),-np.sin(Angle1),0],[np.sin(Angle1),np.cos(Angle1),0],[0,0,1]]) #rotation matrix for rotation about z axis by Theta1 rad
    RAngle4z=np.array([[np.cos(Angle4),-np.sin(Angle4),0],[np.sin(Angle4),np.cos(Angle4),0],[0,0,1]]) #rotation matrix for rotation about z axis by Theta1 rad
    RAngle5y=np.array([[np.cos(Angle5),0,np.sin(Angle5)],[0,1,0],[-np.sin(Angle5),0,np.cos(Angle5)]]) #rotation matrix for rotation about z axis by Theta1 rad
    if resn=='His':
        A2=np.array([((d2*np.cos(Angle2))),(d2*np.sin(Angle2)),0]) #Position of one C that is neighboring Nep in imidazole ring
        A3=np.array([((d2*np.cos(Angle2))),(-d2*np.sin(Angle2)),0]) #Position of other C that is neighboring Nep in imidazole ring
    if resn=='Asp' or resn=='Glu':
        A2=np.array([d2,0,0])
        A3=np.array([(d2+(d2*np.cos(np.pi-Angle2))),(d2*np.sin(np.pi-Angle2)),0])
    A1Rot1=np.matmul(RAngle3x,A1) #Position of A1 after rotating imidazole ring about x axis (shouldn't change)
    A2Rot1=np.matmul(RAngle3x,A2) #Position of A2 after rotating imidazole ring about x axis
    A3Rot1=np.matmul(RAngle3x,A3) #Position of A3 after rotating imidazole ring about x axis
    A1Rot2=np.matmul(RAngle4z,A1Rot1) #Position of A1 after rotating imidazole ring about z axis (shouldn't change)
    A2Rot2=np.matmul(RAngle4z,A2Rot1) #Position of A2 after rotating imidazole ring about z axis
    A3Rot2=np.matmul(RAngle4z,A3Rot1) #Position of A3 after rotating imidazole ring about z axis
    A1Rot3=np.matmul(RAngle5y,A1Rot2) #Position of A1 after rotating imidazole ring about y axis (shouldn't change)
    A2Rot3=np.matmul(RAngle5y,A2Rot2) #Position of A2 after rotating imidazole ring about y axis
    A3Rot3=np.matmul(RAngle5y,A3Rot2) #Position of A3 after rotating imidazole ring about y axis
    A1Trans1=A1Rot3+np.array([d1,0,0]) #Translate to generate metal-ligand bond distance
    A2Trans1=A2Rot3+np.array([d1,0,0]) #Translate to generate metal-ligand bond distance
    A3Trans1=A3Rot3+np.array([d1,0,0]) #Translate to generate metal-ligand bond distance
    A1Rot4=np.matmul(RAngle1z,A1Trans1) #Position of A1 after rotating imidazole ring 109.5-90 rad about z axis
    A2Rot4=np.matmul(RAngle1z,A2Trans1) #Position of A2 after rotating imidazole ring 109.5-90 rad about z axis
    A3Rot4=np.matmul(RAngle1z,A3Trans1) #Position of A3 after rotating imidazole ring 109.5-90 rad about z axis
    return A1Rot4, A2Rot4, A3Rot4

#Orients the protein based on the 3 coordinates obtained in previous function
def orient_chain_from_parameters(inputpose,resi,atom1,atom2,atom3,NepCoord,C1Coord,C2Coord):
    Translation_Ref_Old=inputpose.residue(resi).xyz(atom1) #Vector to move pose such that Nep is at (0,0,0)
    print(Translation_Ref_Old)
    NumpyTransVector_Ref_Old=np.array([-Translation_Ref_Old.x,-Translation_Ref_Old.y,-Translation_Ref_Old.z]) #Converting pyrosetta vector to numpy array
    print(NumpyTransVector_Ref_Old)
    for r in range(1,inputpose.total_residue() + 1):
        for a in range(1,inputpose.residue(r).natoms()+1):
            numpy_coord=np.array([inputpose.residue(r).xyz(a).x,inputpose.residue(r).xyz(a).y,inputpose.residue(r).xyz(a).z]) #Gets the xyz coordinate of atom a of the residue r and converts it to a numpy array
            translate=numpy_coord+NumpyTransVector_Ref_Old
            pyrosetta_translated_coord=xyzVector_double_t() #Now convert the updated coordinate from numpy array to pyrosetta vector
            pyrosetta_translated_coord.x=translate[0]
            pyrosetta_translated_coord.y=translate[1]
            pyrosetta_translated_coord.z=translate[2]
            inputpose.residue(r).set_xyz(a,pyrosetta_translated_coord) #update the coordinate of atom a of residue r.
    #The following three vectors define the reference frame of the protein before alignment. They are used to create a basis set matrix.
    NumpyOldLocalVector_Nep_C1=np.array([(inputpose.residue(resi).xyz(atom1).x+inputpose.residue(resi).xyz(atom2).x),(inputpose.residue(resi).xyz(atom1).y+inputpose.residue(resi).xyz(atom2).y),(inputpose.residue(resi).xyz(atom1).z+inputpose.residue(resi).xyz(atom2).z)])
    NumpyOldLocalVector_Nep_C2=np.array([(inputpose.residue(resi).xyz(atom1).x+inputpose.residue(resi).xyz(atom3).x),(inputpose.residue(resi).xyz(atom1).y+inputpose.residue(resi).xyz(atom3).y),(inputpose.residue(resi).xyz(atom1).z+inputpose.residue(resi).xyz(atom3).z)])
    NumpyOldLocalVector_Cross=np.cross(NumpyOldLocalVector_Nep_C1,NumpyOldLocalVector_Nep_C2)
    NumpyBasisOld=np.array([[NumpyOldLocalVector_Nep_C1[0],NumpyOldLocalVector_Nep_C2[0],NumpyOldLocalVector_Cross[0]],[NumpyOldLocalVector_Nep_C1[1],NumpyOldLocalVector_Nep_C2[1],NumpyOldLocalVector_Cross[1]],[NumpyOldLocalVector_Nep_C1[2],NumpyOldLocalVector_Nep_C2[2],NumpyOldLocalVector_Cross[2]]])
    NumpyBasisOldInv=np.linalg.inv(NumpyBasisOld) #Invert basis set
    #The following three vectors define what the reference frame of the protein should be after alignment (target frame). They are also used to create a basis set matrix
    NewBasisVec1=C1Coord-NepCoord
    NewBasisVec2=C2Coord-NepCoord
    NewBasisVec3=np.cross(NewBasisVec1,NewBasisVec2)
    NumpyBasisNew=np.array([[NewBasisVec1[0],NewBasisVec2[0],NewBasisVec3[0]],[NewBasisVec1[1],NewBasisVec2[1],NewBasisVec3[1]],[NewBasisVec1[2],NewBasisVec2[2],NewBasisVec3[2]]])
    for r in range(1,inputpose.total_residue() + 1):
        for a in range(1,inputpose.residue(r).natoms()+1):
            numpy_coord=np.array([inputpose.residue(r).xyz(a).x,inputpose.residue(r).xyz(a).y,inputpose.residue(r).xyz(a).z]) #Gets the xyz coordinate of atom a of the residue r and converts it to a numpy array
            first_transform=np.matmul(NumpyBasisOldInv,numpy_coord) #applying the inverse of the basis set of the old reference frame converts the atom coordinates from the global frame to the local frame. Note that the local coordinates are the same regardless of which local frame is used.
            second_transform=np.matmul(NumpyBasisNew,first_transform) #applying the basis set of the new reference frame converts the atom coordinates from the local frame back to the global frame, but by using the relationship between the new local frame and the global frame, the coords have been updated to align with the new local frame, which orients the protein correctly
            third_transform=second_transform+NepCoord #Nep is still centered at (0,0,0), so just need to finish by translating everything such that Nep is located at its target position
            pyrosetta_updated_coord=xyzVector_double_t() #Now convert the updated coordinate from numpy array to pyrosetta vector
            pyrosetta_updated_coord.x=third_transform[0]
            pyrosetta_updated_coord.y=third_transform[1]
            pyrosetta_updated_coord.z=third_transform[2]
            inputpose.residue(r).set_xyz(a,pyrosetta_updated_coord) #update the coordinate of atom a of residue r.

#Generates symmetric assembly from single, placed chain
def generate_symmdef_file(Sym,Angle,A1Coord,resi,j,A2,A3,A4):
    Sym_number=int(Sym[1]) #Gets just the integer from the symmetry argument (ex: The 'C3' argument will return 3)
    RThetaSymy=np.array([[np.cos(Angle),0,np.sin(Angle)],[0,1,0],[-np.sin(Angle),0,np.cos(Angle)]]) #rotation matrix for rotation about z axis by ThetaSym rad
    # Store the interfaces in this list. The element order specifies
    # the interface identity (first, second etc) and the value in the
    # list specifies the number of such interfaces. Taken from Rosetta make_symdef_file_denovo.py
    interface_list = []
    if Sym_number == 2:
        interface_list.append(1) #for C2 dimer, there are 2 subunits and 1 unique interface
    if Sym_number == 3:
        interface_list.append(3)
    if Sym_number == 4:
        interface_list.append(4)
        interface_list.append(2) #for C4 tetramer, there are 4 subunits and 2 unique interfaces
    if Sym_number == 5:
        interface_list.append(5)
        interface_list.append(5)
    if Sym_number == 6:
        interface_list.append(6)
        interface_list.append(6)
        interface_list.append(3)
	#	For large rings assume the subunit 1 and subunit >4 do
	#	not interact
    if Sym_number > 6:
        interface_list.append(Sym_number)
        interface_list.append(Sym_number)
        interface_list.append(Sym_number)
    x_vectors=[] #Symdef file needs two vectors for the reference frames of each chain
    y_vectors=[] #These will be contained in two lists of vectors, where element 0 will be x or y vector of reference frame for chain A, element 1 corresponds to chain B, etc.
    x_vectors.append(np.array([1,0,0])) #x vector for ref frame of chain A
    y_vectors.append(np.array([0,1,0])) #y vector for ref frame of chain A
    A1Coord_vectors=[] #Also need position of one atom for each chain - this will be the atom coordinating the metal
    A1Coord_vectors.append(A1Coord)
    for sub in np.arange(0,Sym_number-1):
        x_vectors.append(np.matmul(RThetaSymy,x_vectors[sub])) #x vectors for ref frame of chains B-X
        y_vectors.append(np.matmul(RThetaSymy,y_vectors[sub])) #y vectors for ref frame of chains B-X
        A1Coord_vectors.append(np.matmul(RThetaSymy,A1Coord_vectors[sub])) #Positions of coordinating atoms of chains B-X
    # x3=np.matmul(RThetaSymy,x2) #x vector for ref frame of chain C
    # y3=np.matmul(RThetaSymy,y2) #y vector for ref frame of chain C
    # A1B=np.matmul(RThetaSymy,A1Coord) #position of Nep of chain B
    # A1C=np.matmul(RThetaSymy,A1B) #position of Nep of chain C
    with open('Symdef_Files/Resi'+str(resi)+'_Rotamer'+str(j)+'_Angle3_'+str(A2)+'_Angle4_'+str(A3)+'_Angle5_'+str(A4)+'.sym','w') as w:
        w.write('symmetry_name Resi'+str(resi)+'_Rotamer'+str(j)+'\n')
        # Write the E line using the information stored in interfaces. Also taken from make_symdef_file_denovo.py
        E_string = 'E = '+ str(Sym_number) + '*VRT1'
        interface_num=0
        for interface in interface_list:
            interface_num=interface_num+1
        # disregard interfaces with values == 0
            if interface !=0:
                E_string = E_string + ' + ' + str(interface) + '*(VRT1:VRT'+str(interface_num+1)+')'
        w.write(E_string+'\n')
        w.write('anchor_residue '+str(resi)+'\n')
        w.write('virtual_coordinates_start'+'\n')
        for chain in np.arange(0,len(x_vectors)):
            w.write('xyz VRT'+str(chain+1)+' '+str(x_vectors[chain][0])+','+str(x_vectors[chain][1])+','+str(x_vectors[chain][2])+' '+str(y_vectors[chain][0])+','+str(y_vectors[chain][1])+','+str(y_vectors[chain][2])+' '+str(A1Coord_vectors[chain][0])+','+str(A1Coord_vectors[chain][1])+','+str(A1Coord_vectors[chain][2])+'\n')
        w.write('xyz VRT'+str(Sym_number+1)+' 1.000,0.000,0.000 0.000,1.000,0.000 0.000,0.000,0.000'+'\n')
        w.write('virtual_coordinates_stop'+'\n')
        jump_count=0
        for num_chains in np.arange(0,len(x_vectors)):
            jump_count=jump_count+1
            w.write('connect_virtual JUMP_'+str(jump_count)+' VRT'+str(num_chains+1)+' SUBUNIT'+'\n')
            jump_count=jump_count+1
            w.write('connect_virtual JUMP_'+str(jump_count)+' VRT'+str(Sym_number+1)+' VRT'+str(num_chains+1)+'\n')


def main(args):
    
    import os
    import pandas as pd
    import pyrosetta
    from pyrosetta import init
    from pyrosetta.rosetta.core.import_pose import pose_from_file
    from pyrosetta.rosetta.core.pose import Pose
    from pyrosetta.rosetta.core.scoring import ScoreFunction
    from pyrosetta.rosetta.protocols.loops import get_fa_scorefxn
    from pyrosetta.rosetta.protocols.simple_moves import MutateResidue
    from pyrosetta.rosetta.core.pack.task import TaskFactory
    from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
    from pyrosetta.rosetta.core.pack import rotamer_set, create_packer_graph
    from pyrosetta.rosetta.protocols.simple_moves import SwitchResidueTypeSetMover
    from pyrosetta.rosetta.protocols.symmetry import SetupForSymmetryMover
    from pyrosetta.rosetta.basic import options
    from pyrosetta.rosetta.core.scoring.sasa import SasaCalc
    
    init()
    
    His=pyrosetta.rosetta.core.chemical.AA.aa_his
    Asp=pyrosetta.rosetta.core.chemical.AA.aa_asp
    Glu=pyrosetta.rosetta.core.chemical.AA.aa_glu

    startingpdb=args.pdb
    options.set_real_option('docking:dock_mcm_trans_magnitude',3)
    start=Pose()
    pose_from_file(start, startingpdb)
    BeginningMutateMover=MutateResidue() #Mover that will make a single point mutation
    
    # Make standard score function (ScoreFunction) and a centroid score function (CenScoreFunction)
    ScoreFunction=ScoreFunction()
    ScoreFunction=get_fa_scorefxn()
    CenScoreFunction=pyrosetta.create_score_function("interchain_cen")
    
    #Convert  theta1 and theta_sym to rad
    theta1_rad=(args.theta1-90)*((2*np.pi)/360)
    theta_sym_rad=args.theta_sym*((2*np.pi)/360)
    
    if args.resn=='His':
        resn='His'
        d_internal=1.3
        A_internal=(54)*((2*np.pi)/360)
        atom1,atom2,atom3=10,8,9
        BeginningMutateMover.set_res_name(His) #Mover that will make a single point mutation. This mover is used to make a point mutation to generate the rotamer library for a residue.
        BeginningMutateMover.set_target(1)
        BeginningMutateMover.apply(start) #Mutates resi1 to His. This enables the rotamers of metal binding residue to be obtained later by calling residue 1. If necessary, just mutate residue 1 back to original residue after script is done
        MainMutateMover=MutateResidue() #Mover that will make a single point mutation. This mover is used to make a point mutation to generate the metal binding site.
        MainMutateMover.set_res_name(His)

    if args.resn=='Asp':
        resn='Asp'
        d_internal=1.26
        A_internal=(123)*((2*np.pi)/360)
        atom1,atom2,atom3=8,6,7
        BeginningMutateMover.set_res_name(Asp)
        BeginningMutateMover.set_target(1)
        BeginningMutateMover.apply(start) #Mutates resi1 to Asp
        MainMutateMover=MutateResidue() #Mover that will make a single point mutation. This mover is used to make a point mutation to generate the metal binding site.
        MainMutateMover.set_res_name(Asp)
    
    if args.resn=='Glu':
        resn='Glu'
        d_internal=1.26
        A_internal=(123)*((2*np.pi)/360)
        atom1,atom2,atom3=9,7,8
        BeginningMutateMover.set_res_name(Glu)
        BeginningMutateMover.set_target(1)
        BeginningMutateMover.apply(start) #Mutates resi1 to Glu
        MainMutateMover=MutateResidue() #Mover that will make a single point mutation. This mover is used to make a point mutation to generate the metal binding site.
        MainMutateMover.set_res_name(Glu)


    try:
        os.mkdir('Symdef_Files')
    except:
        pass

    try:
        os.mkdir('Metal_Docker')
    except:
        pass

    #This generates a new pose, sets up and applies the packer to the pose, and get sthe rotamer library from the packer.
    #This pakerpose never gets used again. Its sole purpose is to obtain a rotamer library for a single histidine residue.
    packerpose=Pose()
    packerpose.assign(start)
    TF=TaskFactory()
    PT=TF.create_packer_task(packerpose)
    PT.restrict_to_repacking()
    PT.set_bump_check(False)
    PRM=PackRotamersMover(ScoreFunction,PT)
    PRM.apply(packerpose)
    RSF=rotamer_set.RotamerSetFactory()
    RS=RSF.create_rotamer_set(packerpose)
    RS.set_resid(1)
    PG=create_packer_graph(packerpose,ScoreFunction,PT)
    RS.build_rotamers(packerpose,ScoreFunction,PT,PG) #Gets a library of rotamers for the resi specified by RS.set_resid() (in this case, resi 73, which is a histidine, so it generates a rotamer library for histidine)
    chi_list=[]
    for i in range(1, RS.num_rotamers()+1):
        chi_list.append(RS.rotamer(i).chi()) #Gets the chi angles for each rotamer in the rotamer library and stores them in a list


    mutatepose=Pose()
    orientpose=Pose()
    sympose=Pose()
    cenpose=Pose()
    switch=SwitchResidueTypeSetMover("centroid")
    SasaCalc=SasaCalc()
    pdblist=[]
    SASAlist=[]
    censcorelist=[]
    for resi in args.pos:
        mutatepose.assign(start)
        MainMutateMover.set_target(resi)
        MainMutateMover.apply(mutatepose) #Mutates target resi to His
        for j in range(1, RS.num_rotamers()+1):
            orientpose.assign(mutatepose)
            orientpose.residue(resi).set_all_chi(chi_list[j-1]) #Sets the chi angles of the histidine to those of rotamer j from the rotamer library
            for A2 in args.theta2:
                Angle2rad=A2*((2*np.pi)/360)
                for A3 in args.theta3:
                    Angle3rad=A3*((2*np.pi)/360) 
                    for A4 in args.theta4:
                        Angle4rad=A4*((2*np.pi)/360) 
                        A1Coord,A2Coord,A3Coord=generate_coords_from_parameters(args.d1, d_internal, theta1_rad, A_internal, Angle2rad, Angle3rad, Angle4rad, resn) #Generate the coordinates of the histidine Nepsilon and neighboring carbons from the input parameters
                        orient_chain_from_parameters(orientpose, resi, atom1, atom2, atom3, A1Coord, A2Coord, A3Coord) #Positions one chain based on coordinates calculated for histidine imidazole atoms
                        sympose.assign(orientpose)
                        generate_symmdef_file(args.sym,theta_sym_rad, A1Coord, resi,j,A2,A3,A4)
                        symmetrize = SetupForSymmetryMover('Symdef_Files/Resi'+str(resi)+'_Rotamer'+str(j)+'_Angle3_'+str(A2)+'_Angle4_'+str(A3)+'_Angle5_'+str(A4)+'.sym')
                        symmetrize.apply(sympose) #Generates symmetric assembly based on symmetry definition file
                        SASA=SasaCalc.calculate(sympose)
                        cenpose.assign(sympose) 
                        switch.apply(cenpose)  #Switch pose representation to centroid mode
                        if CenScoreFunction(cenpose)<args.score_cutoff:
                            sympose.dump_pdb('Metal_Docker/Resi'+str(resi)+'_Rotamer'+str(j)+'_Angle3_'+str(A2)+'_Angle4_'+str(A3)+'_Angle5_'+str(A4)+'.pdb') #Only dump pdb if centroid score is below some cutoff (Note, negative score is favorable). May adjust scoring criteria in future
                            pdblist.append('Resi'+str(resi)+'_Rotamer'+str(j)+'_Angle3_'+str(A2)+'_Angle4_'+str(A3)+'_Angle5_'+str(A4))
                            SASAlist.append(SASA)
                            censcorelist.append(CenScoreFunction(cenpose))
    DockingResults=pd.DataFrame({'File Name':pdblist[0:],'SASA':SASAlist[0:],'Centroid Score':censcorelist[0:]})
    DockingResults.to_csv('Metal_Directed_Docking_Results.csv',index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run Metal-Directed Protein Docking')
    parser.add_argument('--pdb', type=str, help='File path of pdb file of monomeric protein building block.')
    parser.add_argument('--resn', type=str, default='His',help='Three-letter code of identity of the metal-binding residue. His=Histidine, Asp=Aspartic Acid, Glu=Glutamic Acid')
    parser.add_argument('--pos', nargs='+', type=int, help='Positions on protein to place metal-binding residue. No limit in number of values within length of protein.')
    parser.add_argument('--sym', type=str, help='Type of symmetry to generate. Currently only supports circular symmetry, may update later. Ex: C2, C3, C4, etc.')
    parser.add_argument('--d1', type=float, help='Metal-ligand bond distance. Only one value accepted')  
    parser.add_argument('--theta1', type=float, help='Ligand-metal-ligand bond angle in degreees. This is the angle from a (real or virtual) axial ligand. Only one value accepted')
    parser.add_argument('--theta2', nargs='+', type=float, help='First rotational degree of freedom in degrees. No limit in number of values.')
    parser.add_argument('--theta3', nargs='+', type=float, help='Second rotational degree of freedom in degrees. No limit in number of values.')
    parser.add_argument('--theta4', nargs='+', type=float, help='Third rotational degree of freedom in degrees. No limit in number of values.')
    parser.add_argument('--theta_sym', type=float, help='Symmetry angle in degreees. Only one value accepted')
    parser.add_argument('--score_cutoff', type=float, default=0, help='Rosetta centroid score cutoff for deposited designs. Default is any designs with Rosetta centroid score <0 REU will be deposited. Ex: -25 would result in only designs <-25 REU being deposited')
    args = parser.parse_args()
    main(args)
