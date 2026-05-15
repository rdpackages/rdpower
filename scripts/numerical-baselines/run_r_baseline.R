#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

script_path <- tryCatch(normalizePath(sys.frame(1)$ofile, mustWork = TRUE), error = function(e) NA_character_)
if (!is.na(script_path)) {
  repo_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)
} else {
  repo_root <- normalizePath(".", mustWork = TRUE)
}

output <- if (length(args) >= 1) args[[1]] else file.path(repo_root, "docs", "audit", "baselines", "r-current.json")
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(rdrobust)
})

source(file.path(repo_root, "R", "rdpower", "R", "rdpower_fun.R"))
source(file.path(repo_root, "R", "rdpower", "R", "rdpower_package.R"))
source(file.path(repo_root, "R", "rdpower", "R", "rdpower.R"))
source(file.path(repo_root, "R", "rdpower", "R", "rdsampsi.R"))
source(file.path(repo_root, "R", "rdpower", "R", "rdmde.R"))

num <- function(x) {
  if (length(x) == 0 || is.null(x) || is.na(x) || !is.finite(x)) return(NA_real_)
  unname(as.numeric(x))
}

quiet <- function(expr) {
  capture.output(value <- force(expr))
  value
}

keep <- function(obj, fields) {
  out <- list()
  for (field in fields) {
    out[[field]] <- num(obj[[field]])
  }
  out
}

summarize_rdpower <- function(obj) {
  keep(obj, c(
    "power.rbc", "se.rbc", "power.conv", "se.conv",
    "sampsi.l", "sampsi.r", "samph.l", "samph.r",
    "N.l", "N.r", "Nh.l", "Nh.r",
    "tau", "bias.l", "bias.r", "Vl.rb", "Vr.rb", "alpha"
  ))
}

summarize_rdsampsi <- function(obj) {
  keep(obj, c(
    "sampsi.h.tot", "sampsi.h.l", "sampsi.h.r", "sampsi.tot",
    "sampsi.h.tot.cl", "sampsi.h.l.cl", "sampsi.h.r.cl", "sampsi.tot.cl",
    "N.l", "N.r", "Nh.l", "Nh.r",
    "bias.l", "bias.r", "var.l", "var.r",
    "samph.l", "samph.r", "tau", "beta", "alpha", "init.cond", "no.iter"
  ))
}

summarize_rdmde <- function(obj) {
  keep(obj, c(
    "mde", "mde.conv", "se.rbc", "se.conv",
    "sampsi.l", "sampsi.r", "samph.l", "samph.r",
    "N.l", "N.r", "Nh.l", "Nh.r",
    "bias.l", "bias.r", "Vl.rb", "Vr.rb", "alpha", "beta"
  ))
}

data <- read.csv(file.path(repo_root, "R", "rdpower_senate.csv"))
z <- data[c("demvoteshfor2", "demmv")]
covs <- data[c("population", "dopen", "dmidterm")]
cluster <- data$state

cases <- list(
  senate_default = list(
    rdpower = summarize_rdpower(quiet(rdpower(data = z, tau = 5))),
    rdsampsi = summarize_rdsampsi(quiet(rdsampsi(data = z, tau = 5))),
    rdmde = summarize_rdmde(quiet(rdmde(data = z)))
  ),
  senate_covs = list(
    rdpower = summarize_rdpower(quiet(rdpower(data = z, tau = 5, covs = covs))),
    rdsampsi = summarize_rdsampsi(quiet(rdsampsi(data = z, tau = 5, covs = covs))),
    rdmde = summarize_rdmde(quiet(rdmde(data = z, covs = covs)))
  ),
  senate_fixed_bandwidth = list(
    rdpower = summarize_rdpower(quiet(rdpower(data = z, tau = 5, h = c(16, 18), b = c(18, 20), all = TRUE))),
    rdsampsi = summarize_rdsampsi(quiet(rdsampsi(data = z, tau = 5, beta = 0.9, samph = c(18, 19), nratio = 0.5, all = TRUE))),
    rdmde = summarize_rdmde(quiet(rdmde(data = z, beta = 0.75, samph = c(12, 13))))
  ),
  senate_cluster = list(
    rdpower = summarize_rdpower(quiet(rdpower(data = z, tau = 5, cluster = cluster))),
    rdsampsi = summarize_rdsampsi(quiet(rdsampsi(data = z, tau = 5, cluster = cluster))),
    rdmde = summarize_rdmde(quiet(rdmde(data = z, cluster = cluster)))
  )
)

result <- list(
  schema_version = 1,
  package = "rdpower",
  language = "r",
  source = "working-tree",
  timestamp_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  environment = list(
    r = R.version.string,
    platform = R.version$platform,
    rdrobust = as.character(packageVersion("rdrobust")),
    rdrobust_source = system.file(package = "rdrobust")
  ),
  cases = cases
)

jsonlite::write_json(result, output, pretty = TRUE, auto_unbox = TRUE, digits = 16, na = "null", null = "null")
cat("Wrote", output, "\n")
