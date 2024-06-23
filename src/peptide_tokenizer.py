from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdChemReactions
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw
IPythonConsole.ipython_useSVG=True
# Example usage
from rdkit.Chem import rdmolops

def convert_aldehydes_to_acids(smiles_list):
    rxn_smarts = '[CX3H1:1](=O)[H].[OH2:2]>>[CX3:1](=O)[O:2]'
    rxn = AllChem.ReactionFromSmarts(rxn_smarts)
    
    acid_smiles_list = []
    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        ps = rxn.RunReactants((mol, Chem.MolFromSmiles('O')))
        if ps:
            product_mol = ps[0][0]  
            Chem.SanitizeMol(product_mol)  
            acid_smiles = Chem.MolToSmiles(Chem.RemoveHs(product_mol), isomericSmiles=True, canonical=True)
            acid_smiles_list.append(acid_smiles)
        else:
            print('Issue with converting C=O into -COOH')
            acid_smiles_list.append(smiles)
    
    return acid_smiles_list

def break_all_NC_O_bonds_and_get_fragments(smiles):
    mol = Chem.MolFromSmiles(smiles)
    pat = Chem.MolFromSmarts('NC=O')
    
    # Finding the largest ring
    def find_largest_ring(mol):
        sssr = Chem.GetSymmSSSR(mol)
        largest_ring = max(sssr, key=len)
        return set(largest_ring)

    largest_ring = find_largest_ring(mol)
    matches = mol.GetSubstructMatches(pat)
    
    emol = Chem.EditableMol(mol)

    bonds_to_break = []

    for match in matches:
        N_idx, C_idx, O_idx = match
        if N_idx in largest_ring and C_idx in largest_ring:
            bonds_to_break.append((N_idx, C_idx))

    # Break bonds 
    for N_idx, C_idx in sorted(bonds_to_break, reverse=True):  # Sort and reverse to avoid indexing issues
        emol.RemoveBond(N_idx, C_idx)


    fragmented_mol = emol.GetMol()
    

    frags = Chem.GetMolFrags(fragmented_mol, asMols=True, sanitizeFrags=True)
    
    fragment_smiles = [Chem.MolToSmiles(frag) for frag in frags]

    return convert_aldehydes_to_acids(fragment_smiles)


def reassemble_fragments(fragments):
    # Example SMARTS for re-forming a peptide bond between two fragments
    rxn_smarts = '[C:1](=O)[O:2].[N:3]>>[C:1](=O)[N:3]'
    rxn = AllChem.ReactionFromSmarts(rxn_smarts)
    
    # Assuming fragments are provided in the correct order and orientation
    mol = Chem.MolFromSmiles(fragments[0])  # Start with the first fragment
    for frag in fragments[1:]:
        frag_mol = Chem.MolFromSmiles(frag)
        # Combine the current molecule with the next fragment
        ps = rxn.RunReactants((mol, frag_mol))
        if ps:
            mol = ps[0][0]  # Update the molecule to the new product
            Chem.SanitizeMol(mol)  # Sanitize the molecule after each reaction step
        else:
            raise ValueError("Reassembly failed at fragment: " + frag)
    
    # Optional: Clean up the molecule's structure
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    AllChem.UFFOptimizeMolecule(mol)
    mol = Chem.RemoveHs(mol)
    
    return Chem.MolToSmiles(mol)


def find_largest_cycle(products):

    mols = []
    for product in products:
        mols.append(product[0])
    largest_cycle_size = 0
    mol_with_largest_cycle = None
    
    for mol in mols:
        sssr = rdmolops.GetSSSR(mol)  # Get the smallest set of smallest rings
        max_ring_size = max((len(ring) for ring in sssr), default=0)
        
        if max_ring_size > largest_cycle_size:
            largest_cycle_size = max_ring_size
            mol_with_largest_cycle = mol
    
    return mol_with_largest_cycle

def aminoacids_to_cyclic_peptide_smiles(fragments):
    original_smiles = reassemble_fragments(fragments)
    intra_rxn = AllChem.ReactionFromSmarts('([C:1][C:2](=[O:6])[O:3].[N:4][C:5])>>[C:1][C:2](=[O:6])[N:4][C:5]')
    aminoacid = Chem.MolFromSmiles(original_smiles)
    products = intra_rxn.RunReactant(aminoacid,0)
    smiles_output = Chem.MolToSmiles(find_largest_cycle(products))
    return smiles_output


if __name__ == '__main__':
    smiles = 'CC(C)C[C@@H]1NC(=O)[C@@H](CC(C)C)NC(=O)[C@H](CCC(=O)N2CCCC2)NC(=O)[C@@H](C)NC(=O)[C@H]2CCCN2C(=O)[C@H](CC(C)C)NC1=O'
    result_fragments = break_all_NC_O_bonds_and_get_fragments(smiles)
    print(result_fragments)
