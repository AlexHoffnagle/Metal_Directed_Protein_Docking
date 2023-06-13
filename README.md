# Metal_Directed_Protein_Docking
Code to run metal-directed protein docking for the de novo design of metalloprotein assemblies with pre-defined protein building blocks. 
See DOI: TBD

Requires Phyton3, numpy, pandas, and Pyrosetta4

The computational design of metalloprotein assemblies as described in DOI: TBD has the following four steps:
  (1) Parameterization of the metal center of interest.
  (2) Metal-directed protein docking
  (3) Interface design
  (4) AlphaFold2 validation

This script accomplishes step 2. It requires a protein building block (input as a pdb file of a monomeric protein); the three-letter code of the metal binding residue of interest (currently supports histidine, aspartate, or glutamate, may update to also support cysteine); positions to attempt to place the metal-binding residue; the symmetry of the system (currently supports any circular symmetry) and the corresponding symmetry angle (i.e. in C4 symmetry, the chains are related by a 90 degree rotation); and 5 geometric parameters that can be used to describe the metal center (metal-ligand bond distance, ligand-metal-ligand bond angle, and the three rotational degrees of freedom of the ligand). These geometric parameters can be obtained either from the structure of an organometallic complex or metalloprotein active site or by evaluating the parameters of an arbitrary metal center through quantum mechanical calculations.

The metal-directed protein docking outputs the pdb files of any docking geometries that satisfy a Rosetta centroid score cutoff (default = 0 REU). It also outputs csv files containing the Rosetta centroid scores and solvent accessible surface areas (SASA) of all docking geometries that satisfy the score cutoff. SASA is inversely related to buried surface area, so a smaller SASA indicates more buried surface area along the resulting protein interfaces. The best docking geometries typically have low (i.e. very negative) Rosetta centroid scores and low SASAs.

To complete the design process, take your favorite docking geometries and perform either Rosetta design or run ProteinMPNN to design the protein interfaces to stablize the assemblies (step 3 above). Then validate the final designs by predicting the structures with AlphaFold2 (step 4). Successful designs should have an all atom RMSD in the metal center of <1.0 Å between the AlphaFold2 prediction and the desired metal center, as well as pLDDT values >90 for most residues.
