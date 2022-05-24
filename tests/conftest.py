from unittest import mock

import pytest
import synapseclient


@pytest.fixture
def genie_config():
    oncotree_link = (
        "http://oncotree.mskcc.org/api/tumorTypes/tree?version=oncotree_2017_06_21"
    )
    config = {
        "vcf2maf": "syn22493903",
        "cna": "syn11600835",
        "bed": "syn11600834",
        "seg": "syn11600836",
        "fusions": "syn11600837",
        "sample": "syn11600838",
        "patient": "syn11600839",
        "patientCounts": "syn11600832",
        "mafSP": "syn11600810",
        "mutationsInCis": "syn11601206",
        "sampleRetraction": "syn11601155",
        "patientRetraction": "syn11600807",
        "vitalStatus": "syn11600808",
        "bedSP": "syn11600809",
        "validationStatus": "syn11601227",
        "errorTracker": "syn11601244",
        "centerMapping": "syn11601248",
        "processTracker": "syn11604890",
        "clinicalSP": "syn11600812",
        "main": "syn7208886",
        "dbMapping": "syn11600968",
        "md": "syn11605077",
        "consortium": "syn11605415",
        "public": "syn11605416",
        "release": "syn11607356",
        "mafinbed": "syn11600814",
        "fileview": "syn11608914",
        "cbs": "syn11600836",
        "vcf": "syn11608914",
        "maf": "syn11608914",
        "centerMaf": "syn12279903",
        "centerMafView": "syn12292501",
        "oncotreeLink": oncotree_link,
        "releaseFolder": "syn17079016",
        "assayinfo": "syn18404286",
        "logs": "syn10155804",
        "center_config": {
            "SAGE": {
                "center": "SAGE",
                "inputSynId": "syn11601335",
                "stagingSynId": "syn11601337",
                "release": True,
                "mutationInCisFilter": "ON",
            },
            "TEST": {
                "center": "TEST",
                "inputSynId": "syn11601340",
                "stagingSynId": "syn11601342",
                "release": True,
                "mutationInCisFilter": "ON",
            },
        },
        "genie_annotation_pkg": "/path/to/nexus",
        "ethnicity_mapping": "syn7434242",
        "race_mapping": "syn7434236",
        "sex_mapping": "syn7434222",
        "sampletype_mapping": "syn7434273",
    }
    return config


@pytest.fixture()
def syn():
    return mock.create_autospec(synapseclient.Synapse)
