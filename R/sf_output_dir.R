## Folder where models are saved
STOCfree_files <- function(out_path = out_path){

  out_path0 <- out_path
  out_path1 <- paste0("./", out_path)

  if(.Platform$OS.type == "windows") out_path <- out_path1
  if(.Platform$OS.type == "unix") out_path <- paste0("/", out_path)

  ## Does the folder already exists?
  ls_dirs <- list.dirs()

  if(!out_path1 %in% ls_dirs){

    dir.create(paste0(getwd(), out_path))

  }

  if(.Platform$OS.type == "windows") return(out_path)
  if(.Platform$OS.type == "unix") return(out_path0)

}




