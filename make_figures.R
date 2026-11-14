#!/usr/bin/env Rscript
# Figures for the MSK-CHORD immortal time manuscript (base R graphics + survival)
# Usage: Rscript make_figures.R  (reads out/analytic_dataset.csv, writes TIFF/PNG/PDF to out/figures/)
suppressPackageStartupMessages(library(survival))
csv <- "out/analytic_dataset.csv"; od <- "out/figures"; dir.create(od, showWarnings = FALSE, recursive = TRUE)
d0 <- read.csv(csv, stringsAsFactors = FALSE)
d  <- subset(d0, included == 1)
d$T <- d$os_months_from_diagnosis; d$L <- d$dx_to_seq_months
d$S <- d$os_months_from_sequencing; d$E <- d$os_status
FONT <- "Liberation Sans"; CUR <- FONT

save_fig <- function(name, w, h, fun) {          # w, h in inches
  CUR <<- FONT
  tiff(file.path(od, paste0(name, ".tif")), width = w, height = h, units = "in",
       res = 600, compression = "lzw", type = "cairo", family = FONT); fun(); dev.off()
  png(file.path(od, paste0(name, ".png")), width = w, height = h, units = "in",
      res = 300, type = "cairo", family = FONT); fun(); dev.off()
  CUR <<- "Helvetica"
  pdf(file.path(od, paste0(name, ".pdf")), width = w, height = h, family = "Helvetica",
      useDingbats = FALSE); fun(); dev.off()
}

km_naive <- function(x) survfit(Surv(T, E) ~ 1, data = x)
km_lt    <- function(x) survfit(Surv(L, T, E) ~ 1, data = x)
km_seq   <- function(x) survfit(Surv(S, E) ~ 1, data = x)
# number at risk at t: entered by t and still under observation at t
nr <- function(x, type, times) sapply(times, function(t) switch(type,
  naive = sum(x$T >= t), lt = sum(x$L <= t & x$T >= t), seq = sum(x$S >= t)))
s_at  <- function(fit, t) summary(fit, times = t, extend = TRUE)$surv

base_axes <- function(xmax, by, xlab, show_xlab = TRUE) {
  axis(1, at = seq(0, xmax, by), cex.axis = 0.85, mgp = c(3, 0.5, 0), tcl = -0.3)
  axis(2, at = seq(0, 1, 0.2), labels = formatC(seq(0, 1, 0.2), format = "f", digits = 1),
       las = 1, cex.axis = 0.85, mgp = c(3, 0.6, 0), tcl = -0.3)
  box(bty = "l")
  if (show_xlab) mtext(xlab, side = 1, line = 1.7, cex = 0.85)
}
# number-at-risk table, written in line units below the x axis
rt <- function(counts, labels, times, top = 3.1, step = 0.9, cex = 0.72, lab_at) {
  mtext("Number at risk", side = 1, line = top, at = lab_at, adj = 1, cex = cex, font = 2)
  for (i in seq_along(counts)) {
    ln <- top + i * step
    mtext(labels[i], side = 1, line = ln, at = lab_at, adj = 1, cex = cex)
    mtext(format(counts[[i]], big.mark = ","), side = 1, line = ln, at = times, cex = cex)
  }
}

# ---------------- Figure 1 ----------------
f1 <- function() {
  layout(matrix(1:2, 1), widths = c(1.05, 1))
  par(mar = c(3.6, 4.2, 1.6, 1), family = CUR, lwd = 0.8, cex = 1)
  x <- pmax(d$L, 0.03)
  br <- 10^seq(log10(0.03), log10(max(x)) + 1e-9, length.out = 41)
  h <- hist(log10(x), breaks = log10(br), plot = FALSE)
  plot(NA, xlim = log10(c(0.03, 600)), ylim = c(0, max(h$counts) * 1.08), axes = FALSE,
       xlab = "", ylab = "")
  rect(h$breaks[-length(h$breaks)], 0, h$breaks[-1], h$counts, col = "grey75", border = "grey25", lwd = 0.5)
  tk <- c(0.03, 0.1, 1, 10, 100); axis(1, at = log10(tk), labels = c("0", "0.1", "1", "10", "100"),
       cex.axis = 0.85, mgp = c(3, 0.5, 0), tcl = -0.3)
  axis(1, at = log10(c(seq(0.03, 0.09, 0.01), seq(0.2, 0.9, 0.1), 2:9, seq(20, 90, 10), 200, 300, 400, 500)),
       labels = FALSE, tcl = -0.15)
  axis(2, las = 1, cex.axis = 0.85, mgp = c(3, 0.6, 0), tcl = -0.3); box(bty = "l")
  q <- quantile(d$L, c(.5, .9))
  abline(v = log10(q[1]), lty = 2); abline(v = log10(q[2]), lty = 3)
  yt <- par("usr")[4]
  text(log10(q[1]) + 0.05, yt * 0.95, sprintf("Median %.1f", q[1]), adj = 0, cex = 0.72)
  text(log10(q[2]) + 0.05, yt * 0.80, sprintf("P90 %.1f", q[2]), adj = 0, cex = 0.72)
  mtext("Months (log scale)", 1, line = 1.9, cex = 0.8)
  mtext("Number of patients", 2, line = 2.9, cex = 0.8)
  mtext("a", side = 3, line = 0.3, adj = -0.18, font = 2, cex = 1)

  par(mar = c(3.6, 7.6, 1.6, 1), cex = 1)
  ty <- c("Non-Small Cell Lung Cancer" = "Non-small cell lung", "Pancreatic Cancer" = "Pancreas",
          "Colorectal Cancer" = "Colorectal", "Breast Cancer" = "Breast", "Prostate Cancer" = "Prostate")
  st <- t(sapply(names(ty), function(k) quantile(d$L[d$cancer_type == k], c(.25, .5, .75))))
  o <- order(st[, 2], decreasing = TRUE); st <- st[o, ]; y <- seq_len(nrow(st))
  plot(NA, xlim = c(0, 72), ylim = c(0.5, nrow(st) + 0.5), axes = FALSE, xlab = "", ylab = "")
  segments(st[, 1], y, st[, 3], y, lwd = 1.2)
  segments(st[, 1], y - 0.12, st[, 1], y + 0.12); segments(st[, 3], y - 0.12, st[, 3], y + 0.12)
  points(st[, 2], y, pch = 19, cex = 0.9)
  axis(1, at = seq(0, 72, 12), cex.axis = 0.75, mgp = c(3, 0.5, 0), tcl = -0.3)
  axis(2, at = y, labels = ty[rownames(st)], las = 1, cex.axis = 0.78, tick = FALSE, mgp = c(3, 0.4, 0))
  box(bty = "l")
  mtext("Months (median and IQR)", 1, line = 1.9, cex = 0.8)
  mtext("b", side = 3, line = 0.3, adj = -0.5, font = 2, cex = 1)
}
save_fig("Figure1", 7.2, 3.0, f1)

# ---------------- Figure 2 ----------------
fn <- km_naive(d); fc <- km_lt(d); fs <- km_seq(d)
f2 <- function() {
  par(mar = c(7.2, 8.4, 0.8, 1.6), family = CUR, lwd = 0.8, cex = 1)
  plot(fn, conf.int = FALSE, mark.time = FALSE, lty = 2, lwd = 1.4, xlim = c(0, 144), ylim = c(0, 1),
       axes = FALSE, xaxs = "i", yaxs = "i")
  lines(fc, conf.int = FALSE, mark.time = FALSE, lty = 1, lwd = 1.4)
  lines(fs, conf.int = FALSE, mark.time = FALSE, lty = 3, lwd = 1.6, col = "grey35")
  base_axes(144, 24, "Months")
  mtext("Overall survival", 2, line = 2.6, cex = 0.85)
  m1 <- summary(fc)$table["median"]; m2 <- summary(fn)$table["median"]
  legend("topright", bty = "n", cex = 0.75, seg.len = 3, lwd = c(1.4, 1.4, 1.6),
         lty = c(2, 1, 3), col = c("black", "black", "grey35"),
         legend = c(sprintf("From diagnosis, naive (median %.1f)", m2),
                    sprintf("From diagnosis, left-truncated (median %.1f)", m1),
                    sprintf("From sequencing, as released (median %.1f)", summary(fs)$table["median"])))
  tt <- seq(0, 144, 24)
  rt(list(nr(d, "naive", tt), nr(d, "lt", tt), nr(d, "seq", tt)),
     c("Naive", "Left-truncated", "From sequencing"), tt, top = 3.0, lab_at = -10)
}
save_fig("Figure2", 6.3, 4.6, f2)

# ---------------- Figure 3 ----------------
f3 <- function() {
  layout(matrix(1:2, 1))
  e <- subset(d, seq_timing_group == "Early"); l <- subset(d, seq_timing_group == "Late")
  for (p in 1:2) {
    par(mar = c(6.2, 5.6, 1.6, 1.2), family = CUR, lwd = 0.8, cex = 1)
    fe <- if (p == 1) km_naive(e) else km_lt(e); fl <- if (p == 1) km_naive(l) else km_lt(l)
    plot(fe, conf.int = FALSE, mark.time = FALSE, lty = 1, lwd = 1.4, xlim = c(0, 144), ylim = c(0, 1),
         axes = FALSE, xaxs = "i", yaxs = "i")
    lines(fl, conf.int = FALSE, mark.time = FALSE, lty = 2, lwd = 1.4)
    base_axes(144, 24, "Months since diagnosis")
    mtext("Overall survival", 2, line = 2.6, cex = 0.8)
    legend("topright", bty = "n", cex = 0.72, lty = c(1, 2), lwd = 1.4, seg.len = 3,
           legend = c(sprintf("Early (median %.1f)", summary(fe)$table["median"]),
                      sprintf("Late (median %.1f)", summary(fl)$table["median"])))
    mtext(if (p == 1) "a  Naive" else "b  Left-truncated", side = 3, line = 0.3, adj = 0, font = 2, cex = 0.85)
    tt <- seq(0, 144, 48); ty_ <- if (p == 1) "naive" else "lt"
    rt(list(nr(e, ty_, tt), nr(l, ty_, tt)), c("Early", "Late"), tt, top = 3.0, lab_at = -14, cex = 0.68)
  }
}
save_fig("Figure3", 7.2, 3.6, f3)

# ---------------- Figure 4 ----------------
ty <- c("Non-Small Cell Lung Cancer" = "Non-small cell lung", "Colorectal Cancer" = "Colorectal",
        "Breast Cancer" = "Breast", "Prostate Cancer" = "Prostate", "Pancreatic Cancer" = "Pancreas")
f4 <- function() {
  layout(matrix(1:6, 2, byrow = TRUE))
  for (k in names(ty)) {
    s <- subset(d, cancer_type == k); a <- km_naive(s); b <- km_lt(s)
    par(mar = c(5.6, 5.4, 1.6, 1.0), family = CUR, lwd = 0.8, cex = 1)
    plot(a, conf.int = FALSE, mark.time = FALSE, lty = 2, lwd = 1.3, xlim = c(0, 144), ylim = c(0, 1),
         axes = FALSE, xaxs = "i", yaxs = "i")
    lines(b, conf.int = FALSE, mark.time = FALSE, lty = 1, lwd = 1.3)
    axis(1, at = seq(0, 144, 24), cex.axis = 0.72, mgp = c(3, 0.4, 0), tcl = -0.25)
    axis(2, at = seq(0, 1, 0.2), labels = formatC(seq(0, 1, 0.2), format = "f", digits = 1),
         las = 1, cex.axis = 0.72, mgp = c(3, 0.5, 0), tcl = -0.25); box(bty = "l")
    mtext("Months since diagnosis", 1, line = 1.5, cex = 0.7)
    mtext("Overall survival", 2, line = 2.2, cex = 0.7)
    mtext(sprintf("%s (n = %s)", ty[k], format(nrow(s), big.mark = ",")), 3, line = 0.3, adj = 0,
          font = 2, cex = 0.7)
    tt <- seq(0, 144, 48)
    rt(list(nr(s, "naive", tt), nr(s, "lt", tt)), c("Naive", "Left-trunc."), tt, top = 2.9,
       step = 0.8, lab_at = -24, cex = 0.58)
  }
  par(mar = c(5.4, 2, 1.6, 0.8), cex = 1); plot.new()
  legend("left", bty = "n", lty = c(2, 1), lwd = 1.3, seg.len = 2.5, cex = 0.62, inset = -0.05,
         legend = c("From diagnosis, naive", "From diagnosis, left-truncated"))
}
save_fig("Figure4", 7.2, 6.0, f4)

# ---------------- Supplementary Figure S1 ----------------
sc <- function(time, entry = NULL) {
  m <- if (is.null(entry)) coxph(Surv(d[[time]], E) ~ I(stage_highest_recorded == "IV") + strata(cancer_type), data = d)
       else coxph(Surv(d[[entry]], d[[time]], E) ~ I(stage_highest_recorded == "IV") + strata(cancer_type), data = d)
  exp(c(coef(m), confint(m)))
}
hr <- rbind(sc("T"), sc("T", "L"), sc("S"))
rownames(hr) <- c("From diagnosis, naive", "From diagnosis, left-truncated", "From sequencing, as released")
print(round(hr, 2))
s1 <- function() {
  par(mar = c(3.4, 14.5, 0.8, 9), family = CUR, lwd = 0.8, cex = 1)
  y <- 3:1
  plot(NA, xlim = log(c(1, 4)), ylim = c(0.5, 3.5), axes = FALSE, xlab = "", ylab = "")
  abline(v = 0, lty = 1, lwd = 0.8)
  segments(log(hr[, 2]), y, log(hr[, 3]), y, lwd = 1.3)
  points(log(hr[, 1]), y, pch = 15, cex = 0.8)
  axis(1, at = log(c(1, 1.5, 2, 3, 4)), labels = c("1", "1.5", "2", "3", "4"), cex.axis = 0.85,
       mgp = c(3, 0.5, 0), tcl = -0.3)
  mtext(rownames(hr), side = 2, at = y, las = 1, line = 0.5, cex = 0.8)
  mtext(sprintf("%.2f (%.2f-%.2f)", hr[, 1], hr[, 2], hr[, 3]), side = 4, at = y, las = 1,
        line = 0.5, cex = 0.8)
  mtext("HR (95% CI)", side = 4, at = 3.5, las = 1, line = 0.5, cex = 0.8, font = 2)
  mtext("Hazard ratio, stage IV vs I-III (log scale)", 1, line = 1.9, cex = 0.8)
}
save_fig("FigureS1", 6.3, 2.1, s1)

# ---------------- Supplementary Figure S2 ----------------
M <- t(sapply(names(ty), function(k) {
  s <- subset(d, cancer_type == k); a <- km_naive(s); b <- km_lt(s)
  100 * (s_at(a, c(24, 60, 120)) - s_at(b, c(24, 60, 120)))
}))
rownames(M) <- ty; colnames(M) <- c("24", "60", "120"); print(round(M, 1))
s2 <- function() {
  par(mar = c(3.4, 9.5, 2.0, 5.5), family = CUR, lwd = 0.8, cex = 1)
  nr <- nrow(M); nc <- ncol(M); brk <- seq(0, 30, 2.5)
  pal <- grey(seq(0.97, 0.35, length.out = length(brk) - 1))
  image(1:nc, 1:nr, t(M[nr:1, ]), breaks = brk, col = pal, axes = FALSE, xlab = "", ylab = "")
  abline(h = 0.5 + 0:nr, v = 0.5 + 0:nc, col = "white", lwd = 2)
  vv <- t(M[nr:1, ])
  for (i in 1:nc) for (j in 1:nr) text(i, j, sprintf("%.0f", vv[i, j]), cex = 0.85,
                                       col = if (vv[i, j] > 18) "white" else "black")
  axis(1, at = 1:nc, labels = c("24", "60", "120"), tick = FALSE, cex.axis = 0.85, mgp = c(3, 0.4, 0))
  axis(2, at = 1:nr, labels = rev(rownames(M)), las = 1, tick = FALSE, cex.axis = 0.85, mgp = c(3, 0.4, 0))
  mtext("Months since diagnosis", 1, line = 1.9, cex = 0.8)
  # legend bar
  usr <- par("usr"); x0 <- usr[2] + 0.25; x1 <- x0 + 0.18
  yb <- seq(usr[3] + 0.4, usr[4] - 0.4, length.out = length(brk))
  rect(x0, yb[-length(yb)], x1, yb[-1], col = pal, border = NA, xpd = NA)
  rect(x0, yb[1], x1, yb[length(yb)], border = "grey30", xpd = NA, lwd = 0.6)
  text(x1 + 0.08, yb[c(1, 5, 9, 13)], brk[c(1, 5, 9, 13)], adj = 0, cex = 0.7, xpd = NA)
  text(x0 - 0.05, usr[4] + 0.08, "Percentage\npoints", adj = c(0, 0), cex = 0.68, xpd = NA)
}
save_fig("FigureS2", 5.2, 3.0, s2)
cat("done\n")
