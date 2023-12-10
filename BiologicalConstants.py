atom_radius = {
    'H': 0.53,
    'C': 0.67,
    'N': 0.56,
    'O': 0.48,
    'P': 0.98,
    'S': 1.02,
    'Na': 1.86,
    'K': 2.27,
    'Cl': 1.75,
    'Mg': 0.65,
    'Ca': 0.99,
}

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
