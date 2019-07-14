# -*- coding:UTF-8 -*-

from bs4 import BeautifulSoup
import requests
import pandas as pd
import lxml as lxml
import re as re
import numpy as np


def getMetabolite(database="HMDB", ID="HMDB0000002"):
    if database == "HMDB":
        url = "http://www.hmdb.ca/metabolites"
    else:
        print("Only support HMDB now.\n")
        return()

    url = url + '/' + ID + ".xml"


    headers={
        'User-Agent':'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/75.0.3770.100 Safari/537.36'
    }

    response = requests.get(url, headers=headers)

    if response.status_code == 404:
        print(ID, "may be not in the database, please check it.\n")
        return()

    # try:
    #     response = requests.get(url, headers=headers)
    # except:
    #     print(ID, "may be not in the database, please check it.\n")

    soup = BeautifulSoup(response.content)

    metabolite = soup.find("metabolite")
    contents = list(metabolite.contents)
    children = list(metabolite.children)

    for child in children:
        if (child == "\n"):
            children.remove(child)

    item_name = list()
    item_value = list()

    for item in children:
        item_name = item_name + [item.name]
        if (len(list(item.children)) == 1):
            item_value = item_value + [item.text]
        else:
            item_value = item_value + [item.children]

    re_analysis_index = list()
    for i in range(0, len(item_value), 1):
        if type(item_value[i]) is type(metabolite.children):
            re_analysis_index = re_analysis_index + [i]

    if (len(re_analysis_index) > 0):
        for idx in re_analysis_index:
            print(idx)
            if item_name[idx] in ["secondary_accessions", "synonyms"]:
                temp_value = children[idx].text
                temp_value = re.split("\n", temp_value)
                temp_value = remove_from_list(temp_value, "")

            if item_name[idx] in ["predicted_properties"]:
                temp_value = children[idx].text
                temp_value = re.split(r"\n\n", temp_value)
                temp_value = remove_from_list(temp_value, "")

                for idx2 in range(0, len(temp_value), 1):
                    temp_value[idx2] = temp_value[idx2].strip()
                    temp_value[idx2] = remove_from_list(re.split("\n", temp_value[idx2]), "")

                temp_value = pd.DataFrame(temp_value)
                temp_value.columns = ["kind", "value", "source"]

            if item_name[idx] in ["taxonomy"]:
                temp_value = getTaxonomy(children[idx])

            if item_name[idx] in ["ontology"]:
                temp_value = getOntology(children[idx])

            if item_name[idx] in ["experimental_properties"]:
                temp_value = getExperimental_properties(children[idx])

            if item_name[idx] in ["spectra"]:
                temp_value = children[idx].find_all("spectrum")

                spectra = list()
                for item in temp_value:
                    item = item.text
                    item = re.split("\n", item)
                    item = remove_from_list(item, "")
                    spectra = spectra + [item]
                spectra = pd.DataFrame(spectra)
                spectra.columns = ["type", "spectrum_id"]
                temp_value = spectra

            if item_name[idx] in ["biological_properties"]:
                temp_value = getBiological_properties(children[idx])

            if item_name[idx] in ["normal_concentrations", "abnormal_concentrations"]:
                temp_value = getConcentrations(object = children[idx], which=item_name[idx])

            if item_name[idx] in ["diseases"]:
                temp_value = getDiseases(children[idx])

            if item_name[idx] in ["pdb_id", "phenol_explorer_compound_id", "knapsack_id", 'bigg_id', 'biocyc_id']:
                temp_value = np.nan

            if item_name[idx] in ["general_references"]:
                temp_value = getGeneral_references(children[idx])

            if item_name[idx] in ["protein_associations"]:
                temp_value = getProtein_associations(children[idx])

            item_value[idx] = temp_value


    result = dict(zip(item_name, item_value))

    return(result)

























def remove_from_list(list, object):
    '''
    :param list: A list object
    :param object: Object you want to remove from the list
    :return: A list object
    '''
    for x in list:
        if (x == object):
            list.remove(object)
    return(list)


def getTaxonomy(object):
    description = object.find_all("description")[0].text
    description = description.strip()
    direct_parent = object.find_all("direct_parent")[0].text
    kingdom = object.find_all("kingdom")[0].text
    super_class = object.find_all("super_class")[0].text
    class_new = object.find_all("class")[0].text
    sub_class = object.find_all("sub_class")[0].text
    molecular_framework = object.find_all("molecular_framework")[0].text
    alternative_parents = object.find_all("alternative_parents")[0].text
    alternative_parents = re.split('\n', alternative_parents)
    alternative_parents = remove_from_list(alternative_parents, "")

    substituents = object.find_all("substituents")[0].text
    substituents = re.split('\n', substituents)
    substituents = remove_from_list(substituents, "")

    external_descriptor = object.find_all("external_descriptor")
    external_descriptor1 = list()
    for item in external_descriptor:
        external_descriptor1 = external_descriptor1 + [item.text]
    external_descriptor = external_descriptor1

    result = [description, direct_parent, kingdom, super_class, class_new, sub_class, molecular_framework, alternative_parents, \
              substituents, external_descriptor
              ]

    result = pd.DataFrame(
        [['description', 'direct_parent', 'kingdom', 'super_class', 'class_new', 'sub_class', 'molecular_framework',
         'alternative_parents',
         'substituents', 'external_descriptor'
         ],
        [description, direct_parent, kingdom, super_class, class_new, sub_class, molecular_framework,
         alternative_parents,
         substituents, external_descriptor
         ]]
    )

    return(result)





def getOntology(object):
    object = object.find_all("root")
    ontology = list()
    for item in object:
        item_term = item.find_all("term")
        term = list()
        for a in item_term:
            term = term + [a.text]
        ontology = ontology + [term]
    return(ontology)


def getExperimental_properties(object):
    object = object.find_all("property")
    experimental_properties = list()
    for item in object:
        item = item.text
        item = item.strip()
        item = re.split("\n", item)
        item = remove_from_list(item,  "")
        if(len(item) == 2):
            item = item + [np.nan]
        experimental_properties = experimental_properties + [item]
    experimental_properties = pd.DataFrame(experimental_properties)

    experimental_properties.columns = ["kind", "value", "source"]
    return(experimental_properties)



def getBiological_properties(object):
    cellular_locations = object.find("cellular_locations").text
    cellular_locations = re.split("\n", cellular_locations)
    cellular_locations = remove_from_list(cellular_locations, "")

    biospecimen_locations = object.find("biospecimen_locations").text
    biospecimen_locations = re.split("\n", biospecimen_locations)
    biospecimen_locations = remove_from_list(biospecimen_locations, "")

    tissue_locations = object.find("tissue_locations").text
    tissue_locations = re.split("\n", tissue_locations)
    tissue_locations = remove_from_list(tissue_locations, "")

    pathways = object.find("pathways").text
    pathways = re.split(r"\n\n", pathways)
    pathways = remove_from_list(pathways, "")
    pathways = remove_from_list(pathways, "\n")

    new_pathways = list()
    for item in pathways:
        item = re.split("\n", item)
        item = remove_from_list(item, "")
        new_pathways = new_pathways + [item]

    biological_properties = [cellular_locations, biospecimen_locations, tissue_locations, new_pathways]
    return(biological_properties)




# def getConcentrations(object):
#     object = object.text
#     object = re.split(r"\n\n", object)
#     object = remove_from_list(object, "")
#
#     result = list()
#     for item in object:
#         item = item.split("\n")
#         item = remove_from_list(item, "")
#         result = result + [item]
#
#     result = pd.DataFrame(result)
#     result.columns = ["biospecimen", "concentration_value", "concentration_units", \
#                       'patient_age', "patient_sex", "patient_information", 'comment']
#
#     return(result)


def getConcentrations(object, which = "normal_concentrations"):

    object = list(object)
    object = remove_from_list(object, "\n")
    if(len(object) == 0):
        return(np.nan)

    result = list()
    for item in object:
        biospecimen = item.find("biospecimen").text
        concentration_value = item.find("concentration_value").text
        concentration_units = item.find("concentration_units").text
        if which == "normal_concentrations":
            subject_age = item.find("subject_age").text
            subject_sex = item.find("subject_sex").text
            subject_condition = item.find("subject_condition").text
        else:
            subject_age = item.find("patient_age").text
            subject_sex = item.find("patient_sex").text
            subject_condition = item.find("patient_information").text

        try:
            comment = item.find("comment")
        except:
            comment = np.nan

        try:
            references = item.find("reference_text").text
        except:
            references = np.nan

        item = [biospecimen, concentration_value, concentration_units, subject_age, subject_sex, subject_condition, comment, references]
        result = result + [item]

    result = pd.DataFrame(result)
    if which == "normal_concentrations":
        result.columns = ["biospecimen", "concentration_value", "concentration_units", \
                          'subject_age', "subject_sex", "subject_condition", 'comment', 'references']
    else:
        result.columns = ["biospecimen", "concentration_value", "concentration_units", \
                          'patient_age', "patient_sex", "patient_condition", 'comment', "references"]

    return(result)



def getDiseases(object):
    object = object.find_all("disease")

    result = list()
    for item in object:
        item = item.text
        item = re.split(r"\n\n\n", item)
        item = [re.sub("\n", "", x) for x in item]
        result = result + [item]

    return(result)



def getGeneral_references(object):
    object = object.find_all("reference")
    result = list()
    for item in object:
        item = item.text
        item = re.split(r"\n", item)
        item = remove_from_list(item, "")
        # item = [re.sub("\n", "", x) for x in item]
        result = result + [item]

    result = pd.DataFrame(result)
    result.columns = ["reference_text", 'pubmed_id']
    return(result)


def getProtein_associations(object):
    object = object.find_all("protein")
    result = list()
    for item in object:
        item = item.text
        item = item.strip()
        item = re.split(r"\n", item)
        item = remove_from_list(item, "")
        # item = [re.sub("\n", "", x) for x in item]
        result = result + [item]

    result = pd.DataFrame(result)
    result.columns = ["protein_accession", 'name', 'uniprot_id', 'gene_name', 'protein_type']
    return(result)






