import requests
import time
from typing import *
from creds import Credentials


class Swiss:
    token = Credentials.swiss_api_token

    @classmethod
    def send_call_swiss_api(cls, sequence: str | List[str], title):
        if isinstance(sequence, str):
            # force sequence to be List[str]
            sequence = [sequence]

        response = requests.post(
            "https://swissmodel.expasy.org/automodel",
            headers={"Authorization": f"Token {cls.token}"},
            json={
                "target_sequences":
                    sequence,
                "project_title": title
            })
        return response

    @classmethod
    def fetch_result_from_swiss_api(cls, response, return_one=True):
        # Obtain the project_id from the response created above
        project_id = response.json()["project_id"]

        # And loop until the project completes
        while True:
            # We wait for some time
            time.sleep(10)

            # Update the status from the server
            response = requests.get(
                f"https://swissmodel.expasy.org/project/{project_id}/models/summary/",
                headers={"Authorization": f"Token {cls.token}"})

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
