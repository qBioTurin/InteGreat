from pyDataverse.api import NativeApi, DataAccessApi
from pyDataverse.models import Datafile
import os


# This class is a controller for Dataverse
# 1. __init__(self, base_url, token): initialize the object using a base_url of dataverse implementation
#   and a dataverse api token (e.g. base_url='https://mirri-dataverse.di.unito.it')
# 2. create_dataverse(self): TODO
# 3. create_dataset(self, dataverse, json): create a dataset from the json in the dataverse parameter
class DataverseController:
    def __init__(self, base_url, token):
        self.base_url = base_url
        self.api = NativeApi(self.base_url, token)
        DataAccessApi(base_url)  # Init a Data Access to the API

    def create_dataverse(self):
        pass

    def create_dataset(self, dataverse, json):
        resp = self.api.get_dataverse(dataverse)

        if resp.json()['status'] != 'OK':
            return 'ERROR dataverse'

        resp = self.api.create_dataset(dataverse, json)

        if resp.json()['status'] != 'OK':
            return resp.json()

        return resp.json()["data"]["persistentId"]  # TODO create a human readable string

    def publish_dataset(self, DOI, release_type='major'):
        return self.api.publish_dataset(DOI, release_type=release_type)

    def add_datafile_to_dataset(self, DOI, results_dir):
        for file in os.listdir(results_dir):
            df = Datafile()
            df_filename = f'{results_dir}/{file}' 
            df.set({"pid": DOI, "filename": df_filename})
            resp = self.api.upload_datafile(DOI, df_filename, df.json())

        #TODO change return

        return resp.json()

