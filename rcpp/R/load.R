# -----------------------------------------------------------------------------
# File:    rcpp/R/load.R
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Bayesian local
#          false discovery rate for sparse count data with application to the
#          discovery of hotspots in protein domains. Ann. Appl. Stat.
#          16(3):1459-1475. Tables 1 to 3, Figures S1 to S3.
# Purpose: Loader of the Rcpp version: compiles the C++ files of rcpp/src/
#          and puts each C++ function in place of the R function of the same
#          name (the functions of R/distributions.R, the Gibbs sampler of
#          R/bayes_fdr.R and rmixture() of R/simulate.R). Everything else
#          stays as R/setup.R loaded it and calls the C++ functions.
# Inputs:  rcpp/src/*.cpp, rcpp/src/*.h. The environment variable
#          LOCALFDR_RCPP_CACHE names the folder of the compiled files; without
#          it they are kept in the temporary folder of the machine, so the
#          files are compiled once and loaded from there afterwards.
# Outputs: None.
# Usage:   source("R/setup.R")
#          source("rcpp/R/load.R")
#          fit <- bayes_local_fdr(x, model = "NP", gamma = 1)
# Source:  R/setup.R (the loader of the R version).
# -----------------------------------------------------------------------------

# The C++ files, one per R file with functions in C++.
RCPP_SOURCES <- c("distributions.cpp", "bayes_fdr.cpp", "simulate.cpp")

#' Folder of the compiled C++ files.
#'
#' @return Path: `LOCALFDR_RCPP_CACHE`, or a folder next to R's temporary
#'   folder named after the R installation.
rcpp_cache_dir <- function() {
  installation <- gsub("[^A-Za-z0-9.]+", "-", sub("/lib/R$", "", R.home()))
  default <- file.path(dirname(tempdir()), paste0("localfdr-rcpp", installation))
  Sys.getenv("LOCALFDR_RCPP_CACHE", default)
}

#' Compile the C++ files and return their functions.
#'
#' The files are compiled with `-ffp-contract=off`, so that the compiler
#' rounds after every multiplication and addition as R does, and are linked
#' to the BLAS and LAPACK library of R.
#'
#' @param root Repository root.
#' @param cache Folder of the compiled files.
#' @return Environment with the C++ functions, named as the R functions.
compile_rcpp <- function(root, cache = rcpp_cache_dir()) {
  dir.create(cache, showWarnings = FALSE, recursive = TRUE)
  Rcpp::registerPlugin("localfdr", function() {
    list(env = list(PKG_CXXFLAGS = "-ffp-contract=off", PKG_LIBS = "$(LAPACK_LIBS) $(BLAS_LIBS)"))
  })
  compiled <- new.env()
  for (file in RCPP_SOURCES) {
    Rcpp::sourceCpp(file.path(root, "rcpp", "src", file), env = compiled, cacheDir = cache)
  }
  compiled
}

#' Put the C++ functions in place of the R functions of the same names.
#'
#' A C++ function takes the arguments of its R function; the default values
#' of the R function that C++ cannot state (a vector, a call) are copied.
#'
#' @param compiled Output of `compile_rcpp()`.
#' @param envir Environment of the R functions.
#' @return Names of the functions replaced, invisibly.
replace_with_rcpp <- function(compiled, envir = globalenv()) {
  force(compiled)
  for (name in ls(compiled)) {
    cpp <- compiled[[name]]
    r <- get(name, envir = envir, mode = "function")
    stopifnot(identical(names(formals(cpp)), names(formals(r))))
    formals(cpp) <- formals(r)
    environment(cpp) <- envir
    assign(name, cpp, envir = envir)
  }
  invisible(ls(compiled))
}

RCPP_FUNCTIONS <- replace_with_rcpp(compile_rcpp(REPO_ROOT))

#' Run a script of scripts/ with the functions now loaded.
#'
#' Evaluates every statement of the R script except its `source("R/setup.R")`,
#' which would put the R functions back. The options, seeds and output files
#' are those of the R script.
#'
#' @param name Name of the script, without folder and extension.
#' @return NULL, invisibly.
run_r_script <- function(name) {
  code <- parse(file.path(REPO_ROOT, "scripts", paste0(name, ".R")))
  loader <- vapply(code, identical, logical(1), quote(source("R/setup.R")))
  stopifnot(sum(loader) == 1)
  for (statement in code[!loader]) eval(statement, globalenv())
  invisible(NULL)
}
