% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/regscan.R
\name{regscan}
\alias{regscan}
\title{regscan : main function for genome-wide region-level analysis}
\usage{
regscan(
  phenocov = NULL,
  pheno,
  REGIONinfo,
  geno_type,
  pheno_type,
  covlist = NULL,
  data = NULL,
  SNPinfo = NULL,
  vcfname = NULL,
  machr2 = NULL,
  qcinput = NULL,
  info_score = NULL,
  multiallelic = FALSE,
  multial_nmaxalleles = 2,
  mafcut = 0.05,
  rcut = 0.99,
  firthreg = FALSE,
  MLCheatmap = FALSE,
  regionlist = NULL,
  alltests = FALSE,
  edgecut = 0.5,
  tol = 1e-16,
  qcmachr2 = NULL,
  covout = FALSE,
  LDpruning = TRUE,
  singleSNPall = FALSE,
  verbose = FALSE,
  debug = FALSE,
  parallel = FALSE
)
}
\arguments{
\item{phenocov}{a dataframe containing the phenotype(s) and  covariates if applicable}

\item{pheno}{the name of the phenotype column in phenocov to be analyzed}

\item{REGIONinfo}{dataframe including region definition, with 4 columns: chr”, “region”, “start.bp”, “end.bp”, as provided for example by the output of “BigLD” function from the gpart R package}

\item{geno_type}{(required): geno_type=“D” if genotypes are in allele dosage format , or genotype format (geno_type=“G”) ; this argument required for SKAT/SKATO tests}

\item{pheno_type}{(required): pheno_type=“C” if phenotype is continuous, or pheno_type=“D” if dichotomous}

\item{regionlist}{list of characters with region names}
}
\description{
regscan : main function for genome-wide region-level analysis
}
