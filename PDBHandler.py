import numpy as np
from Bio import PDB
import os

# used for type hinting
from Bio.PDB.Structure import Structure
from Bio.PDB.Residue import Residue
from Bio.PDB.Atom import Atom


class PDBHandler:
    @staticmethod
    def show_3D_from_pdb(filename, mode=None):
        """expects a .pdb files and renders it in html-ipynb"""
        assert ".pdb" in filename, f"ERROR: pdb_render expected a .pdb file. {filename} was given"
        os.system(f"python show_3D.py {filename} --mode {mode}")

    @staticmethod
    def parse_pdb_content(pdb_file_path):
        # Create a PDB parser
        parser = PDB.PDBParser(QUIET=True)

        # Parse the PDB file
        structure = parser.get_structure('pdb_structure', pdb_file_path)
        print(f"{structure.center_of_mass() = }")

        for model in structure:
            print(f"{model.child_dict = }")

            for chain in model:
                print(f"{chain.get_unpacked_list() = }")

                for residue in chain:
                    print(f"{residue}")
                    print(f"{residue.child_dict = }")

                    for atom in residue:
                        print(f"    {atom.name = }")
                        if hasattr(atom, "charge"):
                            print(f"    {atom.charge = }")
                        print(f"    {atom.coord = }")
                        print(f"    {atom.radius = }")

                        print(f"    {atom.altloc = }")

                        if hasattr(atom, "anisou"):
                            print(f"    {atom.anisou = }")

                        if hasattr(atom, "sigatm"):
                            print(f"    {atom.sigatm = }")

                        if hasattr(atom, "siguij"):
                            print(f"    {atom.siguij = }")

                        print(f"    {atom.bfactor = }")

                        print(f"    {atom.level = }")
                        print(f"    {atom.occupancy = }")
                        print(f"    {atom.parent = }")
                        if hasattr(atom, "vector"):
                            print(f"    {atom.vector = }")

                        print(f"    {atom.is_disordered() = }")
                        print("")
                    print("")

    @staticmethod
    def _compute_residue_atoms_bonds_pos(residue: Residue, bond_pos: list[list[float, float], ...]) \
            -> list[list[float, float], ...]:
        """for a given residue, it returns a list of position of all the atoms that participate in a chemical bond"""
        # todo implement OXT terminal
        # todo implement residue to residue binding (peptide bond)
        # key bind to value
        residue_bond_table = {
            "SER": {
                "N": ["CA"],
                "CA": ["C", "CB"],
                "C": ["O", "OG"]
            },
            "PRO": {
                "N": ["CA"],
                "CA": ["C", "CB"],
                "C": ["O"],
                "CB": ["CG"],
                "CG": ["CD"],
                "CD": ["N"]
            },
            "PHE": {
                "N": ["CA"],
                "CA": ["C", "CB"],
                "C": ["O"],
                "CB": ["CG"],
                "CG": ["CD1", "CD2"],
                "CD1": ["CE1"],
                "CD2": ["CE2"],
                "CE1": ["CZ"],
                "CE2": ["CZ"]
            },
            "THR": {
                "N": ["CA"],
                "CA": ["C", "CB"],
                "C": ["O"],
                "CB": ["OG1", "CG2"]
            },
            "ASP": {
                "N": ["CA"],
                "CA": ["C", "CB"],
                "C": ["O"],
                "CB": ["CG"],
                "CG": ["OD1", "OD2"]

            },
            "GLU": {
                "N": ["CA"],
                "CA": ["C", "CB"],
                "C": ["O"],
                "CB": ["CG"],
                "CG": ["CD"],
                "CD": ["OE1", "OE2"]
            },
            "ILE": {
                "N": ["CA"],
                "CA": ["C", "CB"],
                "C": ["O"],
                "CB": ["CG1", "CG2"],
                "CG1": ["CD1"]
            },
            "ASN": {
                "N": ["CA"],
                "CA": ["C", "CB"],
                "C": ["O"],
                "CB": ["CG"],
                "CG": ["OD1", "ND2"]
            },
            "ALA": {
                "N": ["CA"],
                "CA": ["C", "CB"],
                "C": ["O"]
            },
            "CYS": {
                "N": ["CA"],
                "CA": ["C", "CB"],
                "C": ["O"],
                "CB": ["SG"]
            },
            "LYS": {
                "N": ["CA"],
                "CA": ["C", "CB"],
                "C": ["O"],
                "CB": ["CG"],
                "CG": ["CD"],
                "CD": ["CE"],
                "CE": ["NZ"]
            },
            "TRP": {
                "N": ["CA"],
                "CA": ["C", "CB"],
                "C": ["O"],
                "CB": ["CG"],
                "CG": ["CD1", "CD2"],
                "CD1": ["NE1"],
                "NE1": ["CE2"],
                "CD2": ["CE2", "CE3"],
                "CE2": ["CZ2"],
                "CE3": ["CZ3"],
                "CZ2": ["CH2"],
                "CZ3": ["CH2"]
            },
            "VAL": {
                "N": ["CA"],
                "CA": ["C", "CB"],
                "C": ["O"],
                "CB": ["CG1", "CG2"]
            },
            "MET": {
                "N": ["CA"],
                "CA": ["C", "CB"],
                "C": ["O"],
                "CB": ["CG"],
                "CG": ["SD"],
                "SD": ["CE"]
            },
            "ARG": {
                "N": ["CA"],
                "CA": ["C", "CB"],
                "C": ["O"],
                "CB": ["CG"],
                "CG": ["CD"],
                "CD": ["NE"],
                "NE": ["CZ"],
                "CZ": ["NH1", "NH2"]
            },
            "LEU": {
                "N": ["CA"],
                "CA": ["C", "CB"],
                "C": ["O"],
                "CB": ["CG"],
                "CG": ["CD1", "CD2"]
            },
            "TYR": {
                "N": ["CA"],
                "CA": ["C", "CB"],
                "C": ["O"],
                "CB": ["CG"],
                "CG": ["CD1", "CD2"],
                "CD1": ["CE1"],
                "CD2": ["CE2"],
                "CE1": ["CZ"],
                "CE2": ["CZ"],
                "CZ": ["OH"]
            },
            "GLN": {
                "N": ["CA"],
                "CA": ["C", "CB"],
                "C": ["O"],
                "CB": ["CG"],
                "CG": ["CD"],
                "CD": ["OE1", "NE2"]
            },
            "GLY": {
                "N": ["CA"],
                "CA": ["C"],
                "C": ["O"]
            },
            "HIS": {
                "N": ["CA"],
                "CA": ["C", "CB"],
                "C": ["O"],
                "CB": ["CG"],
                "CG": ["ND1", "CD2"],
                "CD2": ["NE2"],
                "ND1": ["CE1"],
                "CE1": ["NE2"]
            }
        }
        bond_dict = residue_bond_table.get(residue.resname)
        if bond_dict is not None:
            for first_atom_name, list_of_second_atoms in bond_dict.items():
                for second_atom_name in list_of_second_atoms:
                    first_atom = residue.child_dict.get(first_atom_name)
                    second_atom = residue.child_dict.get(second_atom_name)
                    assert first_atom is not None and second_atom is not None, f"{first_atom_name = } and {second_atom_name = } for residue {residue.resname} " \
                                                                               f"residue_bond_table and/or pdb file"
                    first_atom_coord = first_atom.coord
                    second_atom_coord = second_atom.coord

                    bond_pos.append([first_atom_coord, second_atom_coord])
        return bond_pos

    @staticmethod
    def _compute_peptide_bond_pos(residues):
        peptide_bond_atom_pos: list[list[float, float, float], ...] = []
        for idx in range(1, len(residues)):
            prev_aa = residues[idx-1].child_dict
            current_aa = residues[idx].child_dict
            c_atom = prev_aa.get("C")
            n_atom = current_aa.get("N")
            assert c_atom is not None and n_atom is not None, f"{residues[idx-1]} at index {idx-1} does not have a C " \
                                                              f"atom or {residues[idx]} at index {idx} does not have a " \
                                                              f"n atom. Check pdb file"

            peptide_bond_atom_pos.append([c_atom.coord, n_atom.coord])
        return peptide_bond_atom_pos

    @staticmethod
    def compute_atom_pos_in_bonds_from_structure(structure: Structure):
        # structure = PDB.PDBParser().get_structure("pdb_structure", path_to_pdb)
        residues = tuple(structure.get_residues())

        peptide_bond_atom_pos = PDBHandler._compute_peptide_bond_pos(residues)
        atoms_bond_pos: list[list[float, float, float], ...] = []
        for residue in residues:
            atoms_bond_pos = PDBHandler._compute_residue_atoms_bonds_pos(residue, atoms_bond_pos)
            """print(f"{residue.get_resname()}: {residue.child_dict}")
            for atom in residue.get_atoms():
                print(f"    atom {atom.get_id()}")"""
        atoms_bond_pos.extend(peptide_bond_atom_pos)
        return atoms_bond_pos

    @staticmethod
    def get_bulk_atoms_coord_from_structure(structure: Structure) -> tuple[list[str], list[np.array], list[Atom]]:
        atom_name_list = []
        atom_coord_list = []
        atom_obj_list = []
        residues = tuple(structure.get_residues())
        for residue in residues:
            for atom_name, atom_obj in residue.child_dict.items():
                atom_name_list.append(atom_name)
                atom_coord_list.append(atom_obj.coord)
                atom_obj_list.append(atom_obj)
        return atom_name_list, atom_coord_list, atom_obj_list

    @staticmethod
    def render_pdb_file(structure):
        bonds_pos = PDBHandler.compute_atom_pos_in_bonds_from_structure(structure)
        atom_name_list, atom_coord_list, atom_obj_list = PDBHandler.get_bulk_atoms_coord_from_structure(structure)
        # draw cylinder on bonds_pos and spheres on atom_coord_list


if __name__ == "__main__":
    pass
