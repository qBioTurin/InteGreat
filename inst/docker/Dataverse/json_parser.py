import json


def create_value_dict(value, multiple, type_class, type_name):
    return {'value': value, 'multiple': multiple, 'typeClass': type_class, 'typeName': type_name}


# This class provides utility functions to work with Json
# 1. __init__(self, values): initialize the Json object using a dictionary called values
# 2. add_strain(self, values): add to the Json object the strain information, if needed
# 3. strains_metadata(self): return the strains metadata
# 4. citation_metadata(self): return the citation metadata
# 5. parsing_dataset(self): return the json parsified
# 6. save_json(self): return the json object (optionally, save the json file)
class Json:
    def __init__(self, values):
        self.accession_number = None
        self.collection = None
        self.strain_id = None
        self.dataset_title = values['dataset_title']
        self.author_name = values['author_name']
        self.author_affiliation = values['author_affiliation']
        self.dataset_contact_name = values['dataset_contact_name']
        self.dataset_contact_email = values['dataset_contact_email']
        self.ds_description = values['dataset_description']
        #self.subject = values['subject']  # this param should be an array

    def add_strain(self, values):
        self.strain_id = str(values['strainid*'])
        self.collection = values['collection*']
        self.accession_number = values['accessionNumber']

    def strains_metadata(self):
        strain_id = create_value_dict(self.strain_id, False, 'primitive', 'strainid')
        collection = create_value_dict(self.collection, False, 'primitive', 'collection')
        accession_number = create_value_dict(self.accession_number, False, 'primitive', 'accessionnumber')
        return {'fields': [strain_id, collection, accession_number], 'displayName': 'Strain'}

    def citation_metadata(self):
        subject = create_value_dict(["Medicine, Health and Life Sciences"], True, 'controlledVocabulary', 'subject')

        ds_description = {'value': [
            {'dsDescriptionValue': create_value_dict(self.ds_description, False, 'primitive', 'dsDescriptionValue')}],
            'typeClass': 'compound', 'multiple': True, 'typeName': 'dsDescription'}

        dataset_contact_name = create_value_dict(self.dataset_contact_name, False, 'primitive', 'datasetContactName')
        dataset_contact_email = create_value_dict(self.dataset_contact_email, False, 'primitive', 'datasetContactEmail')
        dataset_contact = {'value': [{'datasetContactEmail': dataset_contact_email,
                                      'datasetContactName': dataset_contact_name}],
                           'typeClass': 'compound', 'multiple': True, 'typeName': 'datasetContact'}

        author_name = create_value_dict(self.author_name, False, 'primitive', 'authorName')
        author_affiliation = create_value_dict(self.author_affiliation, False, 'primitive', 'authorAffiliation')
        author = {'value': [{'authorName': author_name, 'authorAffiliation': author_affiliation}],
                  'typeClass': 'compound', 'multiple': True, 'typeName': 'author'}

        title = {'value': self.dataset_title, 'multiple': False, 'typeClass': 'primitive', 'typeName': 'title'}

        return {'fields': [title, author, dataset_contact, ds_description, subject],
                'displayName': 'Citation Metadata'}

    # TODO it is possible to generalize this method in order to add strains only if there are strains parameters
    # TODO "metadataLanguage" is equal "mi" but we need to fix this parameters
    def parsing_dataset(self):
        # metadata_block = {'citation': self.citation_metadata(), 'strains': self.strains_metadata()}
        metadata_block = {'citation': self.citation_metadata()}
        dataset_version = {'metadataBlocks': metadata_block}
        #return {'metadataLanguage': 'en', 'datasetVersion': dataset_version}
        return {'datasetVersion': dataset_version}

    def save_json(self, write=False):
        json_file = json.dumps(self.parsing_dataset(), indent=2)

        if write:
            with open('file.json', 'w') as file:
                file.writelines(json_file)

        return json_file

