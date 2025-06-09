from Bio.PDB import PDBList, PDBParser, PDBIO

# Step 1: Download
pdbl = PDBList()
pdbl.retrieve_pdb_file('1HVR', pdir='.', file_format='pdb')

# Step 2: Parse and clean
parser = PDBParser(QUIET=True)
structure = parser.get_structure('1HVR', 'pdb1hvr.ent')

for model in structure:
    for chain in model:
        residues = list(chain.get_residues())
        for res in residues:
            if res.get_resname() in ['CSO', 'XK2', 'PEG', 'SO4', 'CL']:
                print(f"Removing residue {res.get_resname()} at position {res.id}")
                chain.detach_child(res.id)

# Step 3: Save
io = PDBIO()
io.set_structure(structure)
io.save('1HVR_cleaned.pdb')
