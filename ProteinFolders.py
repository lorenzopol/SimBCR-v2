import requests
import time
import subprocess
import os
from dataclasses import dataclass

from creds import Credentials
import utils as u


@dataclass
class ESMatlas:
    sequence: str

    def _send_call_esmatlas_api(self):
        command = f'curl -X POST --data "{self.sequence}" https://api.esmatlas.com/foldSequence/v1/pdb/'
        pdb_content = subprocess.run(command, shell=True, capture_output=True, text=True)
        pdb_content = pdb_content.stdout.strip()
        return pdb_content

    @staticmethod
    def _convert_api_response_to_pdb_file(pdb_content: str, pdb_output_filepath: str):
        assert u.is_valid_pdb(
            pdb_content), f"\n[ERROR]: ESMatlas did not return a valid pdb file. \n ===== RETURN OUTPUT =====\n{pdb_content}"

        with open(pdb_output_filepath, "w") as cfile:
            cfile.write(pdb_content)

    def fold(self, pdb_output_filepath: str):
        """fold a given sequence and dump raw pdb content to pdb_output_filepath if fold is successful"""
        pdb_content = self._send_call_esmatlas_api()
        self._convert_api_response_to_pdb_file(pdb_content, pdb_output_filepath)

    def fold_and_show_pdb(self, pdb_output_filepath: str,
                          mode=None, cdr3_range=None):
        """extend self.fold and shows its result"""
        self.fold(pdb_output_filepath)
        u.show_3D_from_pdb(pdb_output_filepath, mode, cdr3_range)


@dataclass
class Swiss:
    """deprecated"""
    # todo, fix sequence to field
    sequence: list[str]
    token = Credentials.swiss_api_token

    def send_call_swiss_api(self, title):
        if isinstance(self.sequence, str):
            # force sequence to be List[str]
            sequence = [self.sequence]
        else:
            sequence = self.sequence

        response = requests.post(
            "https://swissmodel.expasy.org/automodel",
            headers={"Authorization": f"Token {self.token}"},
            json={
                "target_sequences":
                    sequence,
                "project_title": title
            })
        return response

    def fetch_result_from_swiss_api(self, response, return_one=True):
        # Obtain the project_id from the response created above
        project_id = response.json()["project_id"]

        # And loop until the project completes
        while True:
            # We wait for some time
            time.sleep(10)

            # Update the status from the server
            response = requests.get(
                f"https://swissmodel.expasy.org/project/{project_id}/models/summary/",
                headers={"Authorization": f"Token {self.token}"})

            # Update the status
            status = response.json()["status"]

            print('Job status is now', status)

            if status in ["COMPLETED", "FAILED"]:
                break

        response_object = response.json()
        models_url = []
        if response_object['status'] == 'COMPLETED':
            for model in response_object['models']:
                models_url.append(model['coordinates_url'])
        return models_url[0] if return_one and len(models_url) > 0 else models_url

    def fold(self, title):
        """fold a sequence and returns the pdb filepath"""
        response = self.send_call_swiss_api(title)
        model_url = self.fetch_result_from_swiss_api(response)

        # download .gz file
        os.system(f"start {model_url}")

        time.sleep(10)

        look_up_directory = Credentials.download_look_up_directory
        output_filename = u.handle_gzip(look_up_directory)
        print(f"{output_filename = }")

        # run jup-notebook for 3D view
        return output_filename

    def fold_and_show_pdb(self, tile: str,
                          mode=None, cdr3_range=None):
        """extend self.fold and shows its result"""
        output_filename = self.fold(tile)
        u.show_3D_from_pdb(output_filename, mode, cdr3_range)
