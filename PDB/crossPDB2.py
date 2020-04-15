from biopandas.pdb import PandasPdb
import argparse
from argparse import ArgumentParser
import time
#Add some comments
def atom_connections():
    """
    CONECT section in PDB - info about atom connections in "entry" column
    HETATM - atoms not connected to primary structure
    """

    conect_data = pdb.df['OTHERS'][pdb.df['OTHERS']['record_name'] == 'CONECT']
    connections =  conect_data['entry'].str.split()
    hetatms = pdb.df['HETATM']['atom_number'].to_list() #ids of all hetatoms

    atom_connect = list() #atoms that connects to other atoms in primary structure
    atom_noconnect = list() #actual hetatoms

    for atoms in connections:    
        if all(int(ids) in hetatms for ids in atoms):
            atom_noconnect.append(atoms)
        else:
            atom_connect.append(atoms)
    return (atom_connect, atom_noconnect)

def hetatoms_ids():
    """
    Find atoms ids in heratms that connect to atoms from primary structure
    """
    atom_connect, atom_nonconnect = atom_connections()
    hetatms = pdb.df['HETATM']['atom_number'].to_list()

    drop_atom_ids = list()
    nodrop_atom_ids = list()

    for atom in atom_connect:
        atom = list(map(int, atom))
        for idx in atom:
            if idx in hetatms:
                drop_atom_ids.append(idx)
            else:
                nodrop_atom_ids.append(idx)
    return list(set(drop_atom_ids))

def cross_linked_atoms():
    atom_ids = hetatoms_ids()
    hetatoms = pdb.df['HETATM']
    hetatoms = hetatoms[~hetatoms['atom_number'].isin(atom_ids)]
    columns = ['element_symbol','atom_number', 'atom_name', 'residue_name', 'x_coord', 'y_coord', 'z_coord']
    hetatoms = hetatoms[columns]
    return hetatoms


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument('--pdb', '-f', required=True, type=str, help='Input PDB format file')
    args = parser.parse_args()
    
    start = time.time()

    pdb = PandasPdb()
    pdb = pdb.read_pdb(args.pdb)
    cross_atoms = cross_linked_atoms()
    cross_atoms.to_csv('cross_atoms.csv')
    print(cross_atoms)

    end = time.time()
    print("Done! {:.2f} sec".format(end-start))

