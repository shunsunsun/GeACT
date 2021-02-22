#!/usr/bin/env r
#
# A simple example to install one or more packages from Bioconductor

## load docopt and remotes (or devtools) from CRAN
suppressMessages({
  library(docopt)               # we need docopt (>= 0.3) as on CRAN
  library(remotes)              # or can use devtools as a fallback
})

if (!requireNamespace("BiocManager", quietly=TRUE)) {
  stop("Please install 'BiocManager' first, for example via 'install2.r BiocManager'.", call. = FALSE)
}

if (!requireNamespace("stringr", quietly=TRUE)) {
  stop("Please install 'stringr' first, for example via 'install2.r stringr'.", call. = FALSE)
}

## configuration for docopt
doc <- "Usage: installBioc.r [-l LIBLOC] [-h] [-x] [-s] [-d DEPS] [-n NCPUS] [-r RETRYs] [--update] [--error] [--timeout TIMEOUT] [-m MIRROR] [--cran-repos CRANREPOS] [--] [PACKAGES ...]

-l --libloc LIBLOC  location in which to install [default: /usr/local/lib/R/site-library]
-d --deps DEPS       install suggested dependencies as well [default: NA]
-n --ncpus NCPUS     number of processes to use for parallel install [default: getOption]
-r --retry RETRYs    number of retry times while installaton error due to network failure  [default: 5]
-t --timeout TIMEOUT timeout in seconds [default: 180]
-e --error           throw error and halt instead of a warning [default: FALSE]
-s --skipinstalled   skip installing already installed packages [default: FALSE]
-u --update          whether attempt to update old packages [default: FALSE]
-m --mirror MIRROR   Mirror for bioconductor to use [default: getOption]
--cran-repos CRANREPOS  Repos for cran, or NULL for file [default: getOption]
-h --help            show this help text
-x --usage           show help and short example usage"

opt <- docopt(doc)			# docopt parsing

if (opt$usage) {
  cat(doc, "\n\n")
  cat("where PACKAGES... can be one or more Bioconductor package names, or local (binary or source)
package files (where extensions .tar.gz, .tgz and .zip are recognised). Optional
arguments understood by BiocManager::install can be passed interspersed in the PACKAGES, though
this requires use of '--'.

Examples:
  installBioc.r -l /tmp/lib S4Vectors               # install into given library")
  
  q("no")
}

# %%%% parsing extra arguments %%%%
opt$retry <- as.integer(opt$retry)
opt$timeout <- as.integer(opt$timeout)

if (opt$deps == "TRUE" || opt$deps == "FALSE") {
  opt$deps <- as.logical(opt$deps)
} else if (opt$deps == "NA") {
  opt$deps <- NA
} else{
  stop("--deps option not supported")
}

if (opt$ncpus == "getOption") {
  opt$ncpus <- getOption("Ncpus", 1L)
} else if (opt$ncpus == "-1") {
  ## parallel comes with R 2.14+
  opt$ncpus <- max(1L, parallel::detectCores())
} else{
  opt$ncpus <- as.integer(opt$ncpus)
}

if (length(opt$cran_repos) == 1) {
    ## docopt results are characters, so if we meant NULL we have to set NULL
    if (opt$cran_repos == "NULL")  {
        opt$cran_repos <- NULL
    } else {
        if (opt$cran_repos == "getOption") {
            ## as littler can now read ~/.littler.r and/or /etc/littler.r we can preset elsewhere
            opt$cran_repos <- getOption("repos")
        }
    }
}

# %%%% install setup %%%%
## ensure installation is stripped
Sys.setenv("_R_SHLIB_STRIP_"="true")

## set timeout
options(timeout=opt$timeout)

if (opt$mirror != "getOption"){
  options(BioC_mirror=opt$mirror)
}

options(repos = opt$cran_repos)

# %%%% helper functions %%%%
try_install <- function(pkgs, lib, deps = NA, Ncpus = 1, max_try = 10, upd = F, 
                        skip_installed = T, error = F){
  library(stringr)
  success <- F
  try_times <- 1
  
  # get package names if in github mode
  pkgs_reduced <- pkgs
  
  while(!success && try_times <= max_try)
  {
    cat(paste0("Try to install: ", try_times, "\n"))
    # remove already installed packages from pkg vector
    pkgs <- if(skip_installed) setdiff(pkgs, colnames(installed.packages())) else pkgs
    cat(pkgs, "\n")
    
    download_warn <- F
    tryCatch(
      withCallingHandlers({
        BiocManager::install(pkgs, update = upd, ask = F, lib = lib, dependencies = deps, 
                             Ncpus = Ncpus)
        if (!download_warn) success <- T
      }, warning = function(warn){
        cat(warn$message, "\n")
        is_download_error <- grepl("download", warn$message) ||
          grepl("Timeout", warn$message) || grepl("receiving data", warn$message)
        is_meta_error <- grepl("connect error", warn$message)
  		if (is_download_error) {
          if (!download_warn) {
            cat("Downloading error, retrying...\n")
            download_warn <<- T
          }
          invokeRestart("muffleWarning")
        } else{
          if (error) stop(warn$message, call. = F)
        }
      }), 
      error = function(err){
        if (error) stop(err$message, call. = T)
      }, 
      finally = {
        try_times <- try_times + 1
      }
    )
  }
  
  if (!success){
    if(error) stop("Max retries exceeded. Installation failed!!!")
  }
}

try_install(opt$PACKAGES, lib = opt$libloc, dep = opt$deps, Ncpus = opt$ncpus, 
            upd = opt$update, skip_installed = opt$skipinstalled, error = opt$error, 
            max_try = opt$retry)

sapply(list.files(path=tempdir(), pattern="^(repos|libloc).*\\.rds$", full.names=TRUE), unlink)
