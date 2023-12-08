import os
import warnings

import Bio.Seq
import pandas as pd
from Bio.PDB import PDBParser

import ProteinFolders as pf
from PDBHandler import PDBHandler as PDBh
from Sequences import Sequences

from creds import Credentials

warnings.filterwarnings("error")


def generate_new_immuneSIM_dataset():
    # todo: instead of manually modifying main.R, add wrapper to modify from entry point
    path = os.path.join(os.getcwd(), 'from_immuneSIM')
    os.system(f"cd {path} & rscript main.R")


def get_variable_sequence_form_immuneSIM(rebuild=False):
    if rebuild:
        generate_new_immuneSIM_dataset()
    path = os.path.join(os.getcwd(), 'from_immuneSIM')
    df = pd.read_csv(os.path.join(path, "hs_igh_sim.csv"))
    sampled_seq = df.sample()

    sampled_seq_raw_data = sampled_seq.to_dict(orient='records')[0]
    dna_sequence = sampled_seq_raw_data["sequence"]
    aa_sequence = sampled_seq_raw_data["sequence_aa"]

    return dna_sequence, aa_sequence, sampled_seq_raw_data


def calculate_cdr3_range(aa_variable_sequence: str, aa_junction: str) -> tuple[int, int]:
    cdr3_start = aa_variable_sequence.index(aa_junction)
    cdr3_end = cdr3_start + len(aa_junction)
    return cdr3_start, cdr3_end


def convert_cdr3range_to_std(cdr3_start, cdr3_end) -> str:
    return f"{cdr3_start}-{cdr3_end}"


def main():
    dna_variable_sequence, aa_variable_sequence, sampled_seq_raw_data = get_variable_sequence_form_immuneSIM(
        rebuild=False)
    aa_junction = sampled_seq_raw_data["junction_aa"]
    cdr3_start, cdr3_end = calculate_cdr3_range(aa_variable_sequence, aa_junction)
    cdr3_range = convert_cdr3range_to_std(cdr3_start, cdr3_end)
    dna_light_chain = "".join([dna_variable_sequence, str(Sequences.DNA_IGG1_CK).lower()])

    aa_light_chain = Bio.Seq.translate(Bio.Seq.transcribe(dna_light_chain))
    path_to_pdb = os.path.join(os.getcwd(), "pdb_files/first_try.pdb")
    folder = pf.ESMatlas(aa_light_chain)
    folder.fold_and_show_pdb(path_to_pdb, "CDR", cdr3_range)


if __name__ == "__main__":
    main()
