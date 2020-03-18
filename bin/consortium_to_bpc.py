"""Most current consortium release cBioPortal files to BPC Synapse Project"""
import synapseclient
import synapseutils as synu


DATA_FOLDER_SYNID = "syn21241322"
# release_synid = "syn21446275"


def remove_gene_panels(syn, file_mapping, remove_seqassays, remove_centers):
    """Removes gene panels that shouldn't be there"""
    for name in file_mapping:
        gene_name = name.replace("data_gene_panel_", "").replace(".txt", "")
        if gene_name in remove_seqassays or gene_name.startswith(tuple(remove_centers)):
            print(name)
            print(file_mapping[name])
            syn.delete(file_mapping[name])


def main(release_synid):
    """Updates BPC project"""
    syn = synapseclient.login()

    # Get existing BPC cBioPortal release files
    existing_files = syn.getChildren(DATA_FOLDER_SYNID)
    existing_map = {exist['name']: exist['id']
                    for exist in existing_files}
    # Get existing gene panels
    genepanel_folder_synid = existing_map['gene_panels']
    existing_gene_panels = syn.getChildren(genepanel_folder_synid)
    genepanel_map = {exist['name']: exist['id']
                    for exist in existing_gene_panels}
    # Get existing case lists
    caselist_folder_synid = existing_map['case_lists']
    case_list = syn.getChildren(caselist_folder_synid)
    caselist_map = {case['name']: case['id'] for case in case_list}
    
    # Get release files
    release_files = syn.getChildren(release_synid)
    synid_map = {release['name']: release['id'] for release in release_files}

    # Copy gene panel files    
    for name in synid_map:
        if name.startswith("data_gene_panel_"):
            ent = syn.get(synid_map[name], followLink=True,
                          downloadFile=False)
            synu.copy(syn, ent, genepanel_folder_synid,
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
        synu.copy(syn, ent, caselist_folder_synid,
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
        exclude = name.startswith(("data_gene_panel_", "data_clinical.txt",
                                   "case_lists")) or name.endswith(".html")
        if not exclude:
            ent = syn.get(synid_map[name], followLink=True,
                          downloadFile=False)
            synu.copy(syn, ent, DATA_FOLDER_SYNID,
                      setProvnance=None,
                      updateExisting=True)
    

if __name__ == "__main__":
    main()