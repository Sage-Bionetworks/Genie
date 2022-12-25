#suppressPackageStartupMessages(library(synapserutils))
suppressPackageStartupMessages(library(argparse))
parser <- ArgumentParser()
parser$add_argument("release",
                    help = "Release version (ie. 5.3-consortium)")
parser$add_argument("--template_path",
                    help = "Path to dashboardTemplate.Rmd",
                    default="dashboardTemplate.Rmd")
parser$add_argument("--staging",
                    action = "store_true",
                    help = "Use staging files")
parser$add_argument("--testing",
                    action = "store_true",
                    help = "Use testing files")
parser$add_argument("--syn_user",
                    help = "Synapse username")
parser$add_argument("--syn_pass",
                    help = "Synapse password")

args <- parser$parse_args()
genie_user <- args$syn_user
genie_pass <- args$syn_pass
release <- args$release
template_path <- args$template_path
suppressPackageStartupMessages(library(synapser))
suppressPackageStartupMessages(library(rmarkdown))

tryCatch({
  synLogin()
}, error = function(err) {
  synLogin(genie_user, genie_pass)
})

if (args$staging) {
 database_synid_mappingid = 'syn12094210'
} else if (args$testing) {
 database_synid_mappingid = 'syn11600968'
} else{
 database_synid_mappingid = 'syn10967259'
}

# For testing only
# database_synid_mappingid = 'syn10967259'
# genie_user=NULL
# genie_pass=NULL
# release = "9.3-consortium"
# template_path = "genie/dashboardTemplate.Rmd"
#


rmarkdown_path = sprintf('%s.Rmd', release)
file.copy(template_path, rmarkdown_path, overwrite = T)
rmarkdown::render(rmarkdown_path,
                  params = list("genieUser" = genie_user,
                                "geniePass" = genie_pass,
                                "database_synid_mappingid" = database_synid_mappingid,
                                "release" = release))
# Obtain release folder
database_synid_mapping = synTableQuery(sprintf('select * from %s',
                                               database_synid_mappingid))
database_synid_mappingdf = as.data.frame(database_synid_mapping)
release_folder_fileview_synid = database_synid_mappingdf$Id[
  database_synid_mappingdf$Database == "releaseFolder"
]
release_folder = synTableQuery(sprintf("select id from %s where name = '%s'",
                                       release_folder_fileview_synid, release))
release_folder_synid = release_folder$asDataFrame()$id
#Store the html file.
html_path = sprintf('%s.html', release)
release_html_ent = synStore(File(html_path,
                                 parentId = release_folder_synid))
synStore(Wiki(markdown = sprintf("${preview?entityId=%s}",
                                 release_html_ent$properties$id),
              owner = release_folder_synid))
