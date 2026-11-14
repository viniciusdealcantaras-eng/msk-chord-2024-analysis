#!/usr/bin/env Rscript
# ============================================================================
# Verificacao independente do manuscrito "Vies de tempo imortal / truncamento
# a esquerda" usando a coorte MSK-CHORD (msk_chord_2024) do cBioPortal.
#
# Baixa os arquivos brutos do datahub oficial, recalcula tudo do zero com o
# pacote 'survival' (validado) e imprime uma tabela ALEGADO vs OBTIDO.
#
# Uso:  Rscript verificar_CHORD.R
# Requisitos: R (>=4.0) e o pacote 'survival' (ja incluso no R base).
# ============================================================================

suppressPackageStartupMessages(library(survival))
options(timeout = 300)
MO <- 30.4375  # dias por mes

# ---------------------------------------------------------------- download ---
base <- "https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/msk_chord_2024"
files <- c(pat="data_clinical_patient.txt", samp="data_clinical_sample.txt",
           dx ="data_timeline_diagnosis.txt", spec="data_timeline_specimen.txt")
dir.create("chord_data", showWarnings = FALSE)
for (f in files) {
  dest <- file.path("chord_data", f)
  if (!file.exists(dest)) {
    cat("baixando", f, "...\n")
    download.file(file.path(base, f), dest, quiet = TRUE)
  }
}

rd_clin <- function(f) read.delim(file.path("chord_data", f), skip = 4, header = TRUE,
                                  sep = "\t", quote = "", comment.char = "",
                                  stringsAsFactors = FALSE, check.names = TRUE)
rd_tl   <- function(f) read.delim(file.path("chord_data", f), header = TRUE,
                                  sep = "\t", quote = "", comment.char = "",
                                  stringsAsFactors = FALSE, check.names = TRUE)

pat  <- rd_clin(files["pat"]);  samp <- rd_clin(files["samp"])
dx   <- rd_tl(files["dx"]);     spec <- rd_tl(files["spec"])
cat(sprintf("\nlidos: %d pacientes | %d amostras | %d dx | %d espec\n",
            nrow(pat), nrow(samp), nrow(dx), nrow(spec)))

# --------------------------------------------------------------- checagens ---
spec$START_DATE <- suppressWarnings(as.numeric(spec$START_DATE))
spec_min <- tapply(spec$START_DATE, spec$PATIENT_ID, min, na.rm = TRUE)
cat(sprintf("[CHECK] 1o especime no dia 0 em %.1f%% dos pacientes\n",
            100 * mean(spec_min == 0, na.rm = TRUE)))

# descricao do OS_MONTHS (linha 2 do arquivo de paciente)
os_desc <- readLines(file.path("chord_data", files["pat"]), n = 2)[2]
os_desc <- strsplit(os_desc, "\t")[[1]]
idx <- which(toupper(gsub("#","",strsplit(readLines(file.path("chord_data",files["pat"]),n=5)[5],"\t")[[1]]))=="OS_MONTHS")
cat("[CHECK] definicao OS_MONTHS:\n   ",
    substr(os_desc[idx], 1, 140), "...\n")

# --------------------------------------------------------------- montagem ----
dx$START_DATE <- suppressWarnings(as.numeric(dx$START_DATE))
dx_min <- tapply(dx$START_DATE, dx$PATIENT_ID, min, na.rm = TRUE)   # diagnostico mais precoce
delay  <- data.frame(PATIENT_ID = names(dx_min),
                     delay_months = as.numeric(-dx_min) / MO, stringsAsFactors = FALSE)

pat$OS_MONTHS <- suppressWarnings(as.numeric(pat$OS_MONTHS))
pat$deceased  <- as.integer(grepl("^1", pat$OS_STATUS))

samp <- samp[order(samp$SAMPLE_ID), ]
ctype <- samp[!duplicated(samp$PATIENT_ID), c("PATIENT_ID", "CANCER_TYPE")]

d <- merge(pat[, c("PATIENT_ID","OS_MONTHS","deceased","OS_STATUS",
                   "STAGE_HIGHEST_RECORDED","CURRENT_AGE_DEID")],
           delay, by = "PATIENT_ID", all.x = TRUE)
d <- merge(d, ctype, by = "PATIENT_ID", all.x = TRUE)
d$CURRENT_AGE_DEID <- suppressWarnings(as.numeric(d$CURRENT_AGE_DEID))

n_total <- nrow(d)
d <- d[!is.na(d$delay_months) & !is.na(d$OS_MONTHS), ]
n_neg <- sum(d$delay_months < 0)
d <- d[d$delay_months >= 0, ]
d$T_dx  <- d$delay_months + d$OS_MONTHS   # OS a partir do diagnostico
d$entry <- d$delay_months                 # tempo de entrada (truncamento a esq.)
cat(sprintf("[CHECK] total=%d  excluidos(neg)=%d  analisados=%d  obitos=%d\n",
            n_total, n_neg, nrow(d), sum(d$deceased)))

# ---------------------------------------------------------------- helpers ----
km_med <- function(fit) unname(summary(fit)$table["median"])
km_S   <- function(fit, t) { s <- summary(fit, times = t, extend = TRUE); s$surv[length(s$surv)] }

# --------------------------------------------------- (a) distrib. do atraso ---
q <- quantile(d$delay_months, c(.25,.5,.75,.90))
res <- list()
res[["a) atraso mediana (m)"]]  <- c(4.2,  round(q[2],1))
res[["a) atraso IQR baixo"]]    <- c(1.4,  round(q[1],1))
res[["a) atraso IQR alto"]]     <- c(24.6, round(q[3],1))
res[["a) atraso P90 (m)"]]      <- c(74.9, round(q[4],1))
res[["a) %% atraso > 24m"]]     <- c(25.3, round(100*mean(d$delay_months>24),1))

# --------------------------------------------------- (b) atraso por tipo ------
med_by <- function(ct) round(median(d$delay_months[d$CANCER_TYPE==ct]),1)
res[["b) NSCLC (m)"]]      <- c(1.8,  med_by("Non-Small Cell Lung Cancer"))
res[["b) Pancreas (m)"]]   <- c(2.8,  med_by("Pancreatic Cancer"))
res[["b) Colorretal (m)"]] <- c(4.3,  med_by("Colorectal Cancer"))
res[["b) Mama (m)"]]       <- c(14.5, med_by("Breast Cancer"))
res[["b) Prostata (m)"]]   <- c(18.7, med_by("Prostate Cancer"))

# --------------------------------------------------- imortal person-time ------
res[["   imortal %% do tempo"]] <- c(43.0, round(100*sum(d$delay_months)/sum(d$T_dx),1))

# --------------------------------------------------- (c) KM principal ---------
fit_naive <- survfit(Surv(T_dx, deceased) ~ 1, data = d)               # ignora truncamento
fit_corr  <- survfit(Surv(entry, T_dx, deceased) ~ 1, data = d)        # truncamento a esquerda
fit_seq   <- survfit(Surv(OS_MONTHS, deceased) ~ 1, data = d)          # a partir do sequenciamento
res[["c) mediana dx INGENUA"]]   <- c(86.0, round(km_med(fit_naive),1))
res[["c) mediana dx CORRIGIDA"]] <- c(46.2, round(km_med(fit_corr),1))
res[["c) mediana sequenciam."]]  <- c(49.5, round(km_med(fit_seq),1))
res[["d) S60 dx INGENUA (%%)"]]   <- c(59.7, round(100*km_S(fit_naive,60),1))
res[["d) S60 dx CORRIGIDA (%%)"]] <- c(42.8, round(100*km_S(fit_corr,60),1))

# --------------------------------------------------- (e) diferenca espuria ----
med_delay <- median(d$delay_months)
d$late <- d$delay_months >= med_delay
mnc <- function(sub, lt=FALSE) {
  f <- if (lt) survfit(Surv(entry, T_dx, deceased) ~ 1, data = sub)
       else     survfit(Surv(T_dx, deceased) ~ 1, data = sub)
  round(km_med(f),1)
}
early <- d[!d$late,]; latedf <- d[d$late,]
res[["e) cedo  INGENUA (m)"]]   <- c(69.5,  mnc(early,FALSE))
res[["e) tarde INGENUA (m)"]]   <- c(100.9, mnc(latedf,FALSE))
res[["e) cedo  CORRIGIDA (m)"]] <- c(66.9,  mnc(early,TRUE))
res[["e) tarde CORRIGIDA (m)"]] <- c(36.4,  mnc(latedf,TRUE))

# --------------------------------------------------- (f) Cox -------------------
cx <- d[d$STAGE_HIGHEST_RECORDED %in% c("Stage 1-3","Stage 4"), ]
cx$stage4 <- as.integer(cx$STAGE_HIGHEST_RECORDED == "Stage 4")
c_naive <- coxph(Surv(T_dx, deceased) ~ stage4 + strata(CANCER_TYPE), data = cx)
c_corr  <- coxph(Surv(entry, T_dx, deceased) ~ stage4 + strata(CANCER_TYPE), data = cx)
c_seq   <- coxph(Surv(OS_MONTHS, deceased) ~ stage4 + strata(CANCER_TYPE), data = cx)
hr <- function(m,v) round(exp(coef(m)[v]),2)
res[["f) HR estadio IV INGENUO"]]   <- c(2.94, hr(c_naive,"stage4"))
res[["f) HR estadio IV CORRIGIDO"]] <- c(2.47, hr(c_corr,"stage4"))
res[["f) HR estadio IV SEQUENC."]]  <- c(2.27, hr(c_seq,"stage4"))

# ---------------------------------------------------------------- relatorio ---
cat("\n=================== ALEGADO  vs  OBTIDO ===================\n")
tab <- do.call(rbind, lapply(names(res), function(k){
  a <- res[[k]][1]; o <- res[[k]][2]
  dif <- ifelse(a==0, NA, round(100*abs(o-a)/abs(a),1))
  st  <- ifelse(is.na(dif), "?", ifelse(dif<=5, "OK", ifelse(dif<=15,"~","XX")))
  data.frame(item=k, alegado=a, obtido=o, dif_pct=dif, ok=st, stringsAsFactors=FALSE)
}))
print(tab, row.names = FALSE)
cat("\nLegenda: OK = dentro de 5% | ~ = 5-15% | XX = >15% de diferenca\n")
cat(sprintf("\nRESUMO: %d de %d itens dentro de 5%% do alegado.\n",
            sum(tab$ok=="OK"), nrow(tab)))
