import os

import Bio.Seq
import pandas as pd
import time
import gzip

import APIsHelper as APIh
from PDBHandler import PDBHandler as PDBh
from Sequences import Sequences

from creds import Credentials


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
    dna_sequence = str(sampled_seq["sequence"].iloc[0])
    aa_sequence = str(sampled_seq["sequence_aa"].iloc[0])
    return dna_sequence, aa_sequence


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


def main():
    dna_variable_sequence, aa_variable_sequence = get_variable_sequence_form_immuneSIM(rebuild=False)
    dna_light_chain = "".join([str(Sequences.DNA_IGG1_CK).lower(), dna_variable_sequence])
    print(f"{dna_light_chain = }")

    aa_light_chain = Bio.Seq.translate(Bio.Seq.transcribe(dna_light_chain))
    print(f"{aa_light_chain = }")
    print(f"Init folding")
    pdb_content = APIh.ESMatlas.send_call_esmatlas_api(aa_light_chain)
    print(f"Folding done")
    print(f"Init writing pdb")
    path_to_pdb = os.path.join(os.getcwd(), "pdb_files/first_try.pdb")
    with open(path_to_pdb, "w") as cfile:
        cfile.write(pdb_content)
    print(f"writing pdb done")
    print(f"Init rendering")
    PDBh.show_3D_from_pdb(path_to_pdb)


if __name__ == "__main__":
    main()

