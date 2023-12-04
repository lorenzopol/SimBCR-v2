import os

import Bio.Seq
import pandas as pd
import time
import gzip
from Bio.PDB import PDBParser
from io import StringIO

import APIsHelper as APIh
from PDBHandler import PDBHandler as PDBh
from Sequences import Sequences

from creds import Credentials
import warnings
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


def handle_gzip(look_up_directory):
    files = [f for f in os.listdir(look_up_directory) if
             os.path.isfile(os.path.join(look_up_directory, f)) and f.endswith(".gz")]
    files_sorted_by_time = sorted(files, key=lambda x: os.path.getmtime(os.path.join(look_up_directory, x)),
                                  reverse=True)
    file = files_sorted_by_time[0]

    # dispatch file to cwd
    input_file_path = os.path.join(look_up_directory, file)
    os.rename(input_file_path, input_file_path.replace(" ", ""))
    input_file_path = input_file_path.replace(" ", "")
    output_filename = ".".join(file.split(r"\ ")[-1].split(".")[:-1])

    # dump content
    with gzip.open(input_file_path, 'rb') as f_in:
        with open(output_filename, 'wb') as f_out:
            # Read and decompress the content from the input .gz file and write it to the output file
            f_out.write(f_in.read())
    return output_filename


def send_protein_sequence_to_swiss(protein_seq: str | list[str], title):
    response = APIh.Swiss.send_call_swiss_api(protein_seq, title)
    model_url = APIh.Swiss.fetch_result_from_swiss_api(response)

    # download .gz file
    os.system(f"start {model_url}")

    time.sleep(10)

    look_up_directory = Credentials.download_look_up_directory
    output_filename = handle_gzip(look_up_directory)

    # run jup-notebook for 3D view
    PDBh.show_3D_from_pdb(output_filename)


def is_valid_pdb(input_string):
    try:
        parser = PDBParser()
        pdb_file = StringIO(input_string)
        _ = parser.get_structure('temp_structure', pdb_file)
        return True  # Valid PDB format
    except Exception as e:
        print(f"Not a valid PDB file: {e}")
        return False  # Not a valid PDB format


def calculate_cdr3_range(aa_variable_sequence: str, aa_junction: str) -> tuple[int, int]:
    cdr3_start = aa_variable_sequence.index(aa_junction)
    cdr3_end = cdr3_start + len(aa_junction)
    return cdr3_start, cdr3_end


def main():
    dna_variable_sequence, aa_variable_sequence, sampled_seq_raw_data = get_variable_sequence_form_immuneSIM(
        rebuild=False)
    aa_junction = sampled_seq_raw_data["junction_aa"]
    cdr3_start, cdr3_end = calculate_cdr3_range(aa_variable_sequence, aa_junction)
    dna_light_chain = "".join([dna_variable_sequence, str(Sequences.DNA_IGG1_CK).lower()])

    aa_light_chain = Bio.Seq.translate(Bio.Seq.transcribe(dna_light_chain))
    pdb_content = APIh.ESMatlas.send_call_esmatlas_api(aa_light_chain)
    assert is_valid_pdb(pdb_content), f"\n[ERROR]: ESMatlas did not return a valid pdb file. \n ===== RETURN OUTPUT =====\n{pdb_content}"

    path_to_pdb = os.path.join(os.getcwd(), "pdb_files/first_try.pdb")
    with open(path_to_pdb, "w") as cfile:
        cfile.write(pdb_content)

    PDBh.show_3D_from_pdb(path_to_pdb, "CDR", f"{cdr3_start}-{cdr3_end}")


if __name__ == "__main__":
    main()
