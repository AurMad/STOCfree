## Folder where models are saved
STOCfree_files <- function(out_path = out_path){

  out_path <- paste0("./", out_path)

  ## Does the folder already exists?
  ls_dirs <- list.dirs()

  if(!out_path %in% ls_dirs){

    dir.create(paste0(getwd(), out_path))

  }

  return(out_path)

}




