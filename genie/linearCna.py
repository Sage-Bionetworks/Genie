# from process_functions import *
# logging.basicConfig(level=logging.INFO)
# logger = logging.getLogger(__name__)

# def validateFilename(filePath, center):
#   pass
    
# def process():
#   temp = syn.chunkedQuery('select id, name from file where parentId == "syn5017766"')
#   cbsFiles = [i['file.id'] for i in temp if "solo.grd" in i['file.name']]
#   linearCNA = pd.DataFrame(columns = ["Hugo_symbol"])
#   for entityId in cbsFiles:
#       sologrdent = syn.get(entityId)
#       newname = sologrdent.name.split(".")[0].replace("_solo","")
#       genelogrratio = pd.DataFrame(columns=["Hugo_symbol",newname])
#       sologrd = pd.read_csv(sologrdent.path,sep="\t")
#       for i in sologrd.index:
#           row = sologrd.iloc[i]
#           if row['Genes.all']!='0':
#               genes = re.findall(".+\((.+)\)",row['Genes.all'])[0]
#               genes = set(genes.split(","))
#               temp = pd.DataFrame(index= range(len(genes)),columns=["Hugo_symbol",newname])
#               temp['Hugo_symbol'] = genes
#               temp[newname] = row['log2(ratio)']
#               genelogrratio = genelogrratio.append(temp)
#       linearCNA = pd.merge(linearCNA, genelogrratio,on='Hugo_symbol',how='outer')
#       linearCNA = linearCNA.drop_duplicates()
#   linearCNA.to_csv("data_linear_CNA_GRCC.txt", sep="\t",index=False)


# def validate(syn, filePathList, pool, oncotree_url, fileType, center):
#   pass