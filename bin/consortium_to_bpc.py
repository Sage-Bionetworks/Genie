"""Most current consortium release cBioPortal files to BPC Synapse Project"""
"""Consortium releases to internal BPC page"""
import argparse

import synapseclient
from synapseclient import Folder
import synapseutils as synu


DATA_FOLDER_SYNID = "syn21574209"
RELEASE_TABLE = "syn16804261"


def find_release(syn, release):
    """Finds the Synapse id of a private consortium release folder"""
    release_synid = syn.tableQuery("select distinct(parentId) from {} where "
                                   "release = '{}'".format(RELEASE_TABLE,
                                                           release))
    releasedf = release_synid.asDataFrame()
    if releasedf.empty:
        raise ValueError("Please specify correct release value")
    return releasedf.iloc[0,0]


def remove_gene_panels(syn, file_mapping, remove_seqassays, remove_centers):
    """Removes gene panels that shouldn't be there"""
    for name in file_mapping:
        gene_name = name.replace("data_gene_panel_", "").replace(".txt", "")
        if gene_name in remove_seqassays or gene_name.startswith(tuple(remove_centers)):
            print(name)
            print(file_mapping[name])
            syn.delete(file_mapping[name])


def main(release):
    """Updates BPC project"""
    if release.endswith("1-consortium"):
        raise ValueError("First consortium release are not released")

    syn = synapseclient.login()
    # Finds the synid of the release
    release_synid = find_release(syn, release)

    # Get existing BPC cBioPortal release files
    bpc_folder_ent = syn.store(Folder(release,
                                      parent=DATA_FOLDER_SYNID))
    caselist_folder_ent = syn.store(Folder("case_lists",
                                             parent=bpc_folder_ent))
    genepanel_folder_ent = syn.store(Folder("gene_panels",
                                              parent=bpc_folder_ent))

    # Get existing gene panels
    existing_gene_panels = syn.getChildren(genepanel_folder_ent)
    genepanel_map = {exist['name']: exist['id']
                    for exist in existing_gene_panels}
    # Get existing case lists
    case_list = syn.getChildren(caselist_folder_ent)
    caselist_map = {case['name']: case['id'] for case in case_list}
    
    # Get release files
    release_files = syn.getChildren(release_synid)
    synid_map = {release['name']: release['id'] for release in release_files}

    # Copy gene panel files    
    for name in synid_map:
        if name.startswith("data_gene_panel_"):
            ent = syn.get(synid_map[name], followLink=True,
                          downloadFile=False)
            synu.copy(syn, ent, genepanel_folder_ent.id,
                      setProvnance=None,
                      updateExisting=True)
    # Remove gene panels
    for name in genepanel_map:
        if name not in synid_map:
            print("Removing: {}({})".format(name, genepanel_map[name]))
            syn.delete(genepanel_map[name])


    new_caselists = syn.getChildren(synid_map['case_lists'])
    new_caselist_map = {case['name']: case['id'] for case in new_caselists}
    # Copy case lists
    for name in new_caselist_map:
        ent = syn.get(new_caselist_map[name], followLink=True,
                      downloadFile=False)
        synu.copy(syn, ent, caselist_folder_ent.id,
                  setProvnance=None,
                  updateExisting=True)

    # Remove case lists
    for name in caselist_map:
        if name not in new_caselist_map:
            print("Removing: {}({})".format(name, caselist_map[name]))
            syn.delete(caselist_map[name])
    
    # Copy rest of the files
    for name in synid_map:
        # Do not copy over files with these patterns
        # exclude = name.startswith(("data_gene_panel_", "data_clinical.txt",
        #                            "case_lists")) or name.endswith(".html")
        if not name.startswith(("data_gene_panel_", "data_clinical.txt",
                                "case_lists")):
            ent = syn.get(synid_map[name], followLink=True,
                          downloadFile=False)
            synu.copy(syn, ent, bpc_folder_ent.id,
                      setProvnance=None,
                      updateExisting=True)
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Consortium to BPC')
    parser.add_argument("release",
                        type=str,
                        metavar="8.2-consortium",
                        help="GENIE release version")
    args = parser.parse_args()
    main(args.release)