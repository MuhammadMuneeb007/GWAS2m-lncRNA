options(warn=1)
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 11) stop("Expected 11 arguments")
stats_file <- args[1]
ld_file <- args[2]
vars_file <- args[3]
outdir <- args[4]
n_gwas <- as.numeric(args[5])
n_ref <- as.numeric(args[6])
L <- as.integer(args[7])
coverage <- as.numeric(args[8])
ld_ridge <- as.numeric(args[9])
estimate_s_max_p <- as.integer(args[10])
locus_id <- args[11]

suppressPackageStartupMessages(library(susieR))
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
logcon <- file(file.path(outdir, "susie_rss.log"), open="wt")
sink(logcon, type="output"); sink(logcon, type="message")
on.exit({sink(type="message"); sink(type="output"); close(logcon)}, add=TRUE)

cat("SuSiE-RSS locus:", locus_id, "\n")
cat("GWAS N:", n_gwas, "\n")
cat("Reference N:", n_ref, "\n")
cat("L:", L, " coverage:", coverage, " ridge:", ld_ridge, "\n")

vars <- scan(vars_file, what="character", quiet=TRUE)
p <- length(vars)
if (p < 2) stop("Fewer than two LD variants")
con <- file(ld_file, "rb")
values <- readBin(con, what="numeric", n=p*p, size=4, endian="little")
close(con)
if (length(values) != p*p) stop("LD matrix length mismatch")
R <- matrix(values, nrow=p, ncol=p, byrow=TRUE)
rm(values); gc()
if (any(!is.finite(R))) stop("LD matrix contains non-finite values")
R <- (R + t(R))/2
diag(R) <- 1
if (is.finite(ld_ridge) && ld_ridge > 0) {
  R <- R * (1-ld_ridge)
  diag(R) <- 1
}

stats <- read.delim(stats_file, header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
if (!all(c("REFERENCE_ID","Z_REF") %in% names(stats))) stop("Matched stats missing REFERENCE_ID/Z_REF")
idx <- match(vars, stats$REFERENCE_ID)
if (any(is.na(idx))) stop("LD variants missing from matched GWAS")
stats <- stats[idx,,drop=FALSE]
if (!all(stats$REFERENCE_ID == vars)) stop("GWAS and LD order mismatch")
z <- as.numeric(stats$Z_REF)
if (any(!is.finite(z))) stop("Non-finite Z scores")

s_est <- NA_real_; s_status <- "NOT_RUN"; s_error <- NA_character_
if (p <= estimate_s_max_p) {
  result <- tryCatch(
    list(ok=TRUE, value=susieR::estimate_s_rss(z=z, R=R, n=n_gwas), error=NA_character_),
    error=function(e) list(ok=FALSE, value=NA_real_, error=conditionMessage(e))
  )
  if (result$ok) { s_est <- result$value; s_status <- "SUCCESS" }
  else { s_status <- "FAILED"; s_error <- result$error }
} else {
  s_status <- paste0("SKIPPED_P_GT_", estimate_s_max_p)
}

formals_rss <- names(formals(susieR::susie_rss))
fit_args <- list(
  z=z, R=R, n=n_gwas, L=min(L,p), coverage=coverage,
  estimate_residual_variance=FALSE, estimate_prior_variance=TRUE,
  refine=FALSE, verbose=TRUE
)
if ("estimate_prior_method" %in% formals_rss) fit_args$estimate_prior_method <- "EM"
if ("max_iter" %in% formals_rss) fit_args$max_iter <- 100
if ("tol" %in% formals_rss) fit_args$tol <- 1e-3
if ("R_finite" %in% formals_rss) fit_args$R_finite <- n_ref
fit <- do.call(susieR::susie_rss, fit_args)
saveRDS(fit, file.path(outdir, "susie_fit.rds"))

pip <- fit$pip
if (length(pip) != nrow(stats)) stop("PIP length mismatch")
variant <- stats
variant$PIP <- as.numeric(pip)
variant$IN_95_CREDIBLE_SET <- FALSE
variant$CREDIBLE_SET <- ""

cs_rows <- list()
if (!is.null(fit$sets) && !is.null(fit$sets$cs)) {
  cs <- fit$sets$cs
  for (i in seq_along(cs)) {
    members <- as.integer(cs[[i]])
    members <- members[members >= 1 & members <= nrow(variant)]
    if (length(members) == 0) next
    cs_name <- names(cs)[i]
    if (is.null(cs_name) || is.na(cs_name) || cs_name == "") cs_name <- paste0("CS", i)
    variant$IN_95_CREDIBLE_SET[members] <- TRUE
    for (j in members) {
      existing <- variant$CREDIBLE_SET[j]
      variant$CREDIBLE_SET[j] <- ifelse(existing == "", cs_name, paste(existing, cs_name, sep=";"))
      cs_rows[[length(cs_rows)+1]] <- data.frame(
        LOCUS_ID=locus_id,
        CREDIBLE_SET=cs_name,
        REFERENCE_ID=variant$REFERENCE_ID[j],
        CHR=variant$CHR[j],
        POS=variant$POS[j],
        PIP=variant$PIP[j],
        stringsAsFactors=FALSE
      )
    }
  }
}
variant <- variant[order(-variant$PIP),,drop=FALSE]
write.table(variant, file.path(outdir,"susie_variants.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
if (length(cs_rows) > 0) cs_df <- do.call(rbind, cs_rows) else cs_df <- data.frame(LOCUS_ID=character(),CREDIBLE_SET=character(),REFERENCE_ID=character(),CHR=integer(),POS=integer(),PIP=numeric())
write.table(cs_df, file.path(outdir,"credible_sets.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

summary <- data.frame(
  LOCUS_ID=locus_id,
  N_VARIANTS=p,
  N_CREDIBLE_SETS=ifelse(is.null(fit$sets$cs),0,length(fit$sets$cs)),
  MAX_PIP=max(variant$PIP, na.rm=TRUE),
  CONVERGED=ifelse(is.null(fit$converged),NA,fit$converged),
  SUSIE_L=min(L,p),
  GWAS_N=n_gwas,
  REFERENCE_N=n_ref,
  COVERAGE=coverage,
  LD_RIDGE=ld_ridge,
  stringsAsFactors=FALSE
)
write.table(summary, file.path(outdir,"susie_summary.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

diagnostics <- data.frame(
  LOCUS_ID=locus_id,
  ESTIMATE_S_RSS=s_est,
  ESTIMATE_S_STATUS=s_status,
  ESTIMATE_S_ERROR=s_error,
  MAX_ABS_Z=max(abs(z)),
  stringsAsFactors=FALSE
)
write.table(diagnostics, file.path(outdir,"susie_diagnostics.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
