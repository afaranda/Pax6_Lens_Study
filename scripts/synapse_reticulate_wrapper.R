##############################################################################
#                                                                            #
#  File: synapse_reticulate_wrapper.R                                        #
#  Purpose: Band-aid until SAGE fixes the synapse R client                   #
#  Created: June 12, 2022                                                    #
#  Author: Adam Faranda                                                      #
#                                                                            #
##############################################################################

library(reticulate)
use_condaenv("synapse-env")
synclient <- import("synapseclient")
synLogin <- function(){
  conn <- synclient$login()
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
synStore <- function(entity, used=NULL, executed=NULL){
  x <- conn$store(
    obj=entity,
    used = used,
    executed = executed
  )
  return(x$properties$id)
}

## Wrappers for Activity
Activity <- function(name, description, used=NULL, executed=NULL){
  return(
    synclient$Activity(
      name = "",
      description = "",
      used=used,
      executed=executed
    )
  )
}

## Wrapper for synSetProvenance
synSetProvenance <- function(entity, activity){
  conn$setProvenance(
    entity = entity,
    activity = activity
  )
}


## Wrapper for File constructor
File <- function(path='.', parent = NULL){
  return(
    synclient$File(
      path = path,
      parent = parent
    )
  )
}

## Wrapper for Folder constructor
Folder <- function(name='new', parent = NULL){
  return(
    synclient$Folder(
      name = name,
      parent = parent
    )
  )
}