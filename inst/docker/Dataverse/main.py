from dataverse_controller import DataverseController
from json_parser import Json
import json
import constant
import sys

if __name__ == "__main__":
    dataverse = 'qbiotest'
    dataverse_controller = DataverseController(constant.base_url, constant.token)
    #results_api = api_caller.rest_get('https://susmirri-mbrcapi.di.unito.it/strains')
    json_file = sys.argv[1]
    results_dir = sys.argv[2]

    with open(json_file, 'r') as file:
        values = json.load(file)

    json = Json(values)
    #json.add_strain(results_api[0])
    json.parsing_dataset()
    #json.save_json(True)
    ds_pid = dataverse_controller.create_dataset(dataverse, json.save_json())
    dataverse_controller.add_datafile_to_dataset(ds_pid, results_dir)
    #dataverse_controller.publish_dataset(ds_pid, 'major')

