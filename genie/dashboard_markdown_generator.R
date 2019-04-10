#suppressPackageStartupMessages(library(synapserutils))
suppressPackageStartupMessages(library(argparse))
parser <- ArgumentParser()
parser$add_argument("release",
                    help = "Release version (ie. 5.3-consortium)")
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
suppressPackageStartupMessages(library(synapser))
suppressPackageStartupMessages(library(rmarkdown))

create_markdown <- function(release, release_folder_synid, template) {
  #This function creates the markdown file for the release dashboard
  #params:
  # release:  Release number
  # release_folder_synid:  Synapse id of the release
  # template:  the dashboard template found in this directory read. 

  final_text = ""
  count = 0
  for (line in template) {
    if (final_text != "") {
      fill_in_text = release
      #This count is so that the synId gets replaced at the right sprintf
      if (count == 1) {
        fill_in_text = release_folder_synid
      }
      if (grepl("%s",line)) {
        final_text = paste(final_text, 
                           sprintf(line, fill_in_text),
                           sep = "\n")
        count = count + 1
      } else{
        final_text = paste(final_text, line, sep = "\n")
      }
    } else {
      #This is so that the start of the Rmd isn't an empty line
      final_text = line
    }
  }
  #The file that gets written out is the release.Rmd
  write(final_text, sprintf("%s.Rmd", release))
  return(sprintf("%s.Rmd", release))
}


tryCatch({
  synLogin()
}, error = function(err) {
  synLogin(genie_user, genie_pass)
})
#options(bitmapType='cairo')

if (args$staging) {
 database_synid_mappingid = 'syn12094210'
} else if (args$testing) {
 database_synid_mappingid = 'syn11600968'
} else{
 database_synid_mappingid = 'syn10967259'
}

#Should use this to create a staging version of the dashboard
database_synid_mapping = synTableQuery(sprintf('select * from %s', database_synid_mappingid))
database_synid_mappingdf = as.data.frame(database_synid_mapping)
release_folder_fileview_synid = database_synid_mappingdf$Id[database_synid_mappingdf$Database == "releaseFolder"]

choose_from_release = synTableQuery(paste(sprintf("select distinct(name) as releases from %s", release_folder_fileview_synid),"where name not like 'Release%' and name <> 'case_lists'"))
releases = as.data.frame(choose_from_release)
if (!any(releases$releases %in% release)) {
  stop(sprintf("Must choose correct release: %s", paste0(releases$releases, collapse=", ")))
}

#release_folder_fileview_synid = "syn17019650"
release_folder = synTableQuery(sprintf("select id from %s where name = '%s'", release_folder_fileview_synid, release))
release_folder_synid = release_folder$asDataFrame()$id
template = readLines("dashboardTemplate.Rmd")

rmarkdown_path = create_markdown(release, release_folder_synid, template)
rmarkdown::render(rmarkdown_path,params = list("genieUser"=genie_user, "geniePass"=genie_pass))
#Store the html file.
release_html_ent = synStore(File(sprintf("%s.html",release),parentId=release_folder_synid))
synStore(Wiki(markdown=sprintf("${preview?entityId=%s}",release_html_ent$properties$id),owner=release_folder_synid))
