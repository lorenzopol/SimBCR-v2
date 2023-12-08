import os
import warnings
import argparse

import Bio.Seq
import pandas as pd

import ProteinFolders as pf
from Sequences import Sequences
from InputParser import InputParser

warnings.filterwarnings("error")


def generate_new_immuneSIM_dataset():
    path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'from_immuneSIM')
    os.system(f"cd {path} & rscript main.R")


def get_variable_sequence_form_immuneSIM(rebuild=False):
    if rebuild:
        generate_new_immuneSIM_dataset()
    path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'from_immuneSIM')
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
    parser = argparse.ArgumentParser()
    parser.add_argument("--number_of_seqs", help="[integer], specify the desired number of sequence. ", default=69,
                        required=False)
    parser.add_argument("--species", help="specify the desired specie [hs, mm]", default="hs", required=False)

    # todo add choices = [ig, tr]
    parser.add_argument("--receptor", help="specify the desired receptor [ig, tr]", default="ig", required=False)

    # todo add choices = [h, l, k, a, b]
    parser.add_argument("--chain", help="specify the desired chain [h, l, k, a, b]", default="h", required=False)
    parser.add_argument("--name_repertoire", help="specify the name of the repertoire. If in doubt, leave as default",
                        default=f"NA", required=False)
    args = parser.parse_args()
    mainR_parser = InputParser(args)
    mainR_parser.parse_and_modify_mainR()

    dna_variable_sequence, aa_variable_sequence, sampled_seq_raw_data = get_variable_sequence_form_immuneSIM(
        rebuild=True)

    aa_junction = sampled_seq_raw_data["junction_aa"]
    cdr3_start, cdr3_end = calculate_cdr3_range(aa_variable_sequence, aa_junction)
    cdr3_range = convert_cdr3range_to_std(cdr3_start, cdr3_end)

    # todo: right now even if a tcr-a is simulated, DNA_IGG1_CH1 always gets picked
    dna_light_chain = "".join([dna_variable_sequence, str(Sequences.DNA_IGG1_CH1).lower()])
    aa_light_chain = Bio.Seq.translate(Bio.Seq.transcribe(dna_light_chain))

    path_to_pdb = os.path.join(os.path.dirname(os.path.realpath(__file__)), "pdb_files/first_try.pdb")
    folder = pf.ESMatlas(aa_light_chain)
    folder.fold_and_show_pdb(path_to_pdb, "CDR", cdr3_range)


if __name__ == "__main__":
    main()
