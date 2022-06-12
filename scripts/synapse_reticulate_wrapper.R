##############################################################################
#
#  File: synapse_reticulate_wrapper.R
#  Purpose: Band-aid until SAGE fixes the synapse R client
#  Created: 
#

library(reticulate)
use_condaenv("synapse-env")
sc <- import("synapseclient")
synLogin <- function(){
  conn <- sc$login()
  return(conn)
}
conn <- synLogin()

## Wrapper for synFindEntityId
synFindEntityId <- function(name, parent=NULL){
  return(conn$findEntityId(name = name, parent = parent))
}

## Wrapper for synGet
synGet <- function(entity, downloadLocation='.'){
  conn$get(entity, downloadLocation = downloadLocation)
}

## Wrapper for synStore
synStore <- function(entity){
  x <- conn$store(entity)
  return(x$properties$id)
}

## Wrappers for Activity
Activity <- function(name, description){
  return(
    sc$Activity(
      name = "",
      description = ""
    )
  )
}

syn_act <- Activity(
  name="upload_analysis_results",
  description="upload analysis results"
)
syn_act$executed(syn_script)
## Wrapper for synSetProvenance
synSetProvenance <- function(entity, activity){
  conn$setProvenance(
    entity = entity,
    activity = activity
  )
}

synSetProvenance(x, syn_act)

## Wrapper for File
File <- function(path='.', parent = NULL){
  return(
    sc$File(
      path = path,
      parent = parent
    )
  )
}

