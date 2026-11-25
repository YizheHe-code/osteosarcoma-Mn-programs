# =============================================================================
# 20_supplementary_figures_and_spatial_revision.R
#
# Supplementary Figures S1-S4, the corrected lesion-level spatial analysis and
# the final version of Figure 5. Run from the project root after
# scripts/00_main_all_lessons.R.
#
# Sections
#   1. Shared settings (figure theme taken from 19_revision_checks.R)
#   2. Fig. S1  - single-cell QC and marker-based annotation (discovery cohort)
#   3. Fig. S2  - MR scatter and leave-one-out plots
#   4. Fig. S3  - annotation and tier scores in the validation cohort
#   5. Fig. S4  - spot-level maps of all 26 Visium sections
#   6. Spatial  - lesion-level sensitivity analysis and mixed models (25 lesions)
#   7. Fig. 5   - final figure with panel D based on 25 lesions
#   8. Additional file 1 - merge S1-S5 into one PDF
#
# Note on lesions (sections 6-7): in GSE293065, MR10A and MR10B are two lung
# relapses from the same patient one year apart, whereas MR11A and MR11B are two
# sections of the same metastasis (Moquin-Beaudry et al., Nat Commun 2026;17:8294).
# Lesions are therefore defined by combining MR11A/B only, giving 25 lesions from
# 19 patients.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(data.table); library(ggplot2); library(patchwork)
  library(ggrastr); library(nlme); library(magick); library(qpdf)
})
dir.create("final", showWarnings = FALSE)

# ---- 1. Shared settings ------------------------------------------------------
main = readLines("scripts/00_main_all_lessons.R", encoding = "UTF-8", warn = FALSE)
s19  = readLines("scripts/19_revision_checks.R",  encoding = "UTF-8", warn = FALSE)

# Figure theme and save function (19_revision_checks.R, lines 19-44)
stopifnot(grepl("^save_jtm = function", s19[24]), grepl("^theme_jtm = function", s19[31]))
eval(parse(text = paste(s19[19:44], collapse = "\n")))

save_both = function(p, stem, w, h) {
  ggsave(paste0(stem, ".pdf"),  p, width = w, height = h, units = "mm", device = cairo_pdf)
  ggsave(paste0(stem, ".tiff"), p, width = w, height = h, units = "mm", dpi = 300,
         device = ragg::agg_tiff, compression = "lzw")
}

# ---- 2. Fig. S1: single-cell QC and annotation (GSE152048) --------------------
# Requires `mk` (33 annotation markers; data.table with columns cell_type, gene),
# defined in the cell-type annotation section of 00_main_all_lessons.R.
raw = readRDS("data_proc/sc_raw.rds")
raw[["percent.mt"]] = PercentageFeatureSet(raw, pattern = "^MT-")
md = as.data.table(raw@meta.data, keep.rownames = "cell")[, .(cell, patient, nFeature_RNA, percent.mt)]
rm(raw); invisible(gc())

# QC rule as in the main script: >=300 genes, <= patient median + 3 MAD of log10 genes, <20% mito
md[, lg := log10(nFeature_RNA)]
thr = md[, .(upper = 10^(median(lg) + 3 * mad(lg))), by = patient]
md  = merge(md, thr, by = "patient")
md[, pass := nFeature_RNA >= 300 & nFeature_RNA <= upper & percent.mt < 20]

ann  = readRDS("data_proc/sc_annotated.rds")
fin  = as.data.table(ann@meta.data, keep.rownames = "cell")
drop = c("Lung_epithelial", "Myocyte")                        # 730 cells excluded before scoring
chk  = md[fin, on = "cell"]
stopifnot(sum(!chk$pass) == 0, fin[!celltype %in% drop, .N] == 108863)

lv = paste0("BC", sort(as.integer(sub("BC", "", unique(md$patient)))))
md[, patient := factor(patient, levels = lv)]; thr[, patient := factor(patient, levels = lv)]
xr = theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7))

pA = ggplot(md, aes(patient, nFeature_RNA)) +
  geom_violin(scale = "width", fill = "grey85", colour = "grey40", linewidth = 0.3) +
  geom_hline(yintercept = 300, linetype = 2, colour = "#d1495b") +
  geom_errorbar(data = thr, aes(x = patient, ymin = upper, ymax = upper), inherit.aes = FALSE,
                width = 0.8, colour = "#d1495b") +
  scale_y_log10() + labs(x = NULL, y = "Genes detected per cell") + theme_jtm() + xr
pB = ggplot(md, aes(patient, percent.mt)) +
  geom_violin(scale = "width", fill = "grey85", colour = "grey40", linewidth = 0.3) +
  geom_hline(yintercept = 20, linetype = 2, colour = "#d1495b") +
  labs(x = NULL, y = "Mitochondrial reads (%)") + theme_jtm() + xr
cnt = rbind(md[, .N, by = patient][, stage := "Before QC"],
            fin[!celltype %in% drop, .N, by = patient][, stage := "Analysed"])
cnt[, `:=`(patient = factor(patient, levels = lv), stage = factor(stage, levels = c("Before QC", "Analysed")))]
pC = ggplot(cnt, aes(patient, N, fill = stage)) + geom_col(position = "dodge", width = 0.75) +
  scale_fill_manual(values = c("grey70", "#1b6ca8")) +
  labs(x = NULL, y = "Cells", fill = NULL) + theme_jtm() + theme(legend.position = "top") + xr

genes = unique(mk$gene)
ord = c("OS_Osteoblastic", "OS_Chondroblastic", "Osteoclast", "Myeloid", "T_NK",
        "Fibroblast_MSC", "Endothelial", "Pericyte")
lab = c(OS_Osteoblastic = "OS osteoblastic", OS_Chondroblastic = "OS chondroblastic",
        Osteoclast = "Osteoclast", Myeloid = "Myeloid", T_NK = "T/NK",
        Fibroblast_MSC = "Fibroblast/MSC", Endothelial = "Endothelial", Pericyte = "Pericyte")
dot_theme = theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7, face = "italic"),
                  axis.text.y = element_text(size = 8), legend.text = element_text(size = 7),
                  legend.title = element_text(size = 8))
Idents(ann) = factor(ann$celltype, levels = c(rev(ord), drop))
pD = DotPlot(ann, features = genes, idents = ord, cols = c("grey90", "#1b6ca8")) +
  scale_y_discrete(labels = lab) + labs(x = NULL, y = NULL) + dot_theme   # IGHG1 is absent from the matrix
save_both((pA | pB) / pC / pD + plot_layout(heights = c(1, 0.8, 1.3)) + plot_annotation(tag_levels = "A"),
          "final/FigureS1", 180, 240)
rm(ann, md, chk, fin); invisible(gc())

# ---- 3. Fig. S2: MR scatter and leave-one-out plots ---------------------------
# Requires the FinnGen R13 MR objects from the MR section of the main analysis:
#   iv   instruments with beta_Mn, se_Mn, beta_outcome, se_outcome (one row per SNP x outcome)
#   res  IVW / MR-Egger / weighted-median estimates (columns outcome, method, b)
#   sens Egger intercepts (columns outcome, Egger_intercept)
#   loo  leave-one-out estimates (columns outcome, SNP_removed, OR, OR_lower95, OR_upper95)
short = c("Osteosarcoma (C3_OSTEOSARCOMA_WIDE)" = "Osteosarcoma",
          "Bone and articular cartilage (C3_BONE_CARTILAGE)" = "Bone and articular cartilage")
sc = copy(iv)[, s := sign(beta_Mn)][, `:=`(x = beta_Mn * s, y = beta_outcome * s)]   # Mn-raising allele as effect allele
sc[, out := factor(short[outcome], levels = short)]
ln = merge(res[, .(outcome, method, b)], sens[, .(outcome, Egger_intercept)], by = "outcome")
ln[, a := ifelse(grepl("Egger", method), Egger_intercept, 0)]
ln[, out := factor(short[outcome], levels = short)]
pA = ggplot(sc, aes(x, y)) +
  geom_hline(yintercept = 0, colour = "grey85") +
  geom_segment(aes(x = x, xend = x, y = y - 1.96 * se_outcome, yend = y + 1.96 * se_outcome), colour = "grey65") +
  geom_segment(aes(x = x - 1.96 * se_Mn, xend = x + 1.96 * se_Mn, y = y, yend = y), colour = "grey65") +
  geom_point(size = 1.6) +
  geom_abline(data = ln, aes(intercept = a, slope = b, colour = method, linetype = method), linewidth = 0.6) +
  scale_colour_manual(values = c("#1b6ca8", "#d1495b", "#66a182")) +
  facet_wrap(~out, scales = "free_y") +
  labs(x = "SNP effect on whole-blood Mn (SD)", y = "SNP effect on outcome (log OR)", colour = NULL, linetype = NULL) +
  theme_jtm() + theme(legend.position = "bottom")
lo = copy(loo)[, out := factor(short[outcome], levels = short)]
lo[, all_snp := !grepl("^rs", SNP_removed)]
lo[all_snp == TRUE, SNP_removed := "All SNPs (IVW)"]
lo[, lab := factor(SNP_removed, levels = rev(c(sort(unique(SNP_removed[!all_snp])), "All SNPs (IVW)")))]
pB = ggplot(lo, aes(OR, lab)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_segment(aes(x = OR_lower95, xend = OR_upper95, y = lab, yend = lab, colour = all_snp)) +
  geom_point(aes(colour = all_snp, shape = all_snp), size = 2) +
  scale_colour_manual(values = c(`FALSE` = "grey30", `TRUE` = "#d1495b"), guide = "none") +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 18), guide = "none") +
  scale_x_log10() + facet_wrap(~out, scales = "free_y") +
  labs(x = "OR per SD of whole-blood Mn (95% CI), leaving out the SNP shown", y = NULL) + theme_jtm()
save_both(pA / pB + plot_layout(heights = c(1, 1.1)) + plot_annotation(tag_levels = "A"),
          "final/FigureS2", 180, 210)

# ---- 4. Fig. S3: validation cohort (GSE162454) -------------------------------
vs = readRDS("data_proc/sc_val_annotated.rds")
s1 = as.data.table(openxlsx::read.xlsx("results/tables/Additional_file_2.xlsx", sheet = "S1_tier_genes"))
# Tier gene sets exactly as in the main script (lines 1716-1722); STING1 is TMEM173 in this matrix
stopifnot(grepl("^sting_name", main[1716]), grepl("^mn_sets = lapply", main[1722]))
eval(parse(text = paste(main[1716:1722], collapse = "\n")))
stopifnot(length(setdiff(s1$gene, sub("^TMEM173$", "STING1", unlist(mn_sets)))) == 0)
vs  = AddModuleScore(vs, features = mn_sets, name = "Tier")          # as main script line 1726
tc  = paste0("Tier", 1:4)
vmd = as.data.table(vs@meta.data)[, c("patient", "celltype", tc), with = FALSE]

# Check: patients in whom each cell type exceeds osteoblastic OS (>=20 cells) must match Fig. 3B
pb = vmd[, c(.(n = .N), lapply(.SD, mean)), by = .(patient, celltype), .SDcols = tc][n >= 20]
cc = merge(pb[celltype != "OS_Osteoblastic"], pb[celltype == "OS_Osteoblastic", c("patient", tc), with = FALSE],
           by = "patient", suffixes = c("", "_os"))
print(rbindlist(lapply(tc, function(t) cc[, .(tier = t, k = sum(get(t) > get(paste0(t, "_os"))), N = .N), by = celltype])))

ordv = c("OS_Osteoblastic", "OS_Chondroblastic", "Osteoclast", "Myeloid", "T_NK", "B_Plasma", "Mast",
         "Fibroblast_MSC", "Endothelial", "Pericyte")
labv = c(lab, B_Plasma = "B/plasma", Mast = "Mast")
Idents(vs) = factor(vs$celltype, levels = rev(ordv))
pA = DotPlot(vs, features = unique(c(genes, "TPSAB1", "TPSB2", "CPA3", "KIT")), cols = c("grey90", "#1b6ca8")) +
  scale_y_discrete(labels = labv) + labs(x = NULL, y = NULL) + dot_theme
tl = c(Tier1 = "T1 Mn transport", Tier2 = "T2 Mn enzymes", Tier3 = "T3 Mn glycosylation", Tier4 = "T4 cGAS\u2013STING")
hm = melt(vmd[, lapply(.SD, mean), by = celltype, .SDcols = tc], id.vars = "celltype")
hm[, z := as.numeric(scale(value)), by = variable]
hm[, `:=`(celltype = factor(celltype, levels = rev(ordv)), variable = factor(tl[as.character(variable)], levels = tl))]
pB = ggplot(hm, aes(variable, celltype, fill = z)) + geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.1f", z)), size = 2.5) +
  scale_fill_gradient2(low = "#1b6ca8", mid = "white", high = "#d1495b", name = "z-score") +
  scale_y_discrete(labels = labv) + labs(x = NULL, y = NULL) + theme_jtm() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_both(pA / pB + plot_layout(heights = c(1.4, 1)) + plot_annotation(tag_levels = "A"),
          "final/FigureS3", 180, 200)
rm(vs, vmd); invisible(gc())

# ---- 5. Fig. S4: all 26 Visium sections (GSE293065) --------------------------
# Gene sets (main script lines 1796-1805); scoring and smoothing as main script lines 1903-1937
stopifnot(grepl("^cell_sets *= *list\\(", main[1796]), grepl("^mn_sets *= *list\\(", main[1802]))
eval(parse(text = paste(main[1796:1805], collapse = "\n")))
sp_dir = "data_raw/spatial"
sp_pre = basename(list.dirs(file.path(sp_dir, "samples"), recursive = FALSE))
stopifnot(length(sp_pre) == 26)
sec = sub("^GSM[0-9]+_", "", sp_pre)
grp = ifelse(grepl("^PT", sec), "Primary", ifelse(grepl("^LR", sec), "Local relapse", "Metastasis"))
key = order(factor(grp, levels = c("Primary", "Metastasis", "Local relapse")), as.integer(gsub("\\D", "", sec)), sec)
sp_pre = sp_pre[key]; sec = sec[key]; grp = grp[key]
feat = c("Tumour", "Myeloid", "T2_Enzyme", "Osteoclast", "aVb3")
flab = c("Osteoblastic OS", "Myeloid", "T2 Mn enzymes", "Osteoclast", "\u03b1v\u03b23 integrin")

get_st = function(i) {
  st = NormalizeData(Load10X_Spatial(file.path(sp_dir, "samples", sp_pre[i])), verbose = FALSE)
  show_sets = list(Myeloid = cell_sets$Myeloid, T2_Enzyme = mn_sets$T2_Enzyme,
                   Osteoclast = cell_sets$Osteoclast, Tumour = cell_sets$OS_Osteoblastic)
  show_sets = lapply(show_sets, intersect, rownames(st))
  st = AddModuleScore(st, features = show_sets, name = "Show")
  names(st@meta.data)[grep("^Show", names(st@meta.data))] = names(show_sets)
  ex = GetAssayData(st, layer = "data")
  st$aVb3 = colMeans(as.matrix(ex[intersect(c("ITGAV", "ITGB3"), rownames(ex)), , drop = FALSE]))
  coords = GetTissueCoordinates(st)[colnames(st), ]
  D  = as.matrix(dist(coords[, c("x", "y")]))
  nn = t(apply(D, 1, function(r) order(r)[1:7])); rm(D)        # each spot and its six nearest neighbours
  smooth = function(v) rowMeans(matrix(v[nn], nrow = nrow(nn)))  # display only
  for (f in feat) st[[paste0(f, "_s")]] = smooth(st[[f]][, 1])
  st
}
map_one = function(i) {
  st = get_st(i)
  p  = SpatialFeaturePlot(st, features = paste0(feat, "_s"), ncol = 5, pt.size.factor = 2.5, image.alpha = 0.3)
  for (k in 1:5) p[[k]] = p[[k]] + ggtitle(flab[k]) +
    guides(fill = guide_colourbar(title = NULL, barheight = grid::unit(2.5, "mm"), barwidth = grid::unit(22, "mm"))) +
    theme(plot.title = element_text(size = 8, hjust = 0.5), legend.position = "bottom",
          legend.text = element_text(size = 6), legend.margin = margin(0, 0, 0, 0))
  n_sp[sec[i]] <<- ncol(st)
  lab = wrap_elements(full = grid::textGrob(paste0(sec[i], " (", grp[i], ")"), x = grid::unit(0.01, "npc"),
                                            hjust = 0, gp = grid::gpar(fontsize = 9, fontface = "bold")))
  lab / wrap_elements(full = p) + plot_layout(heights = c(0.07, 1))
}
n_sp = integer(0); pngs = character(0)
pages = split(seq_along(sp_pre), ceiling(seq_along(sp_pre) / 4))
for (j in seq_along(pages)) {
  pl = lapply(pages[[j]], map_one)
  while (length(pl) < 4) pl = c(pl, list(plot_spacer()))
  f = sprintf("final/FigureS4_p%d.png", j)
  ragg::agg_png(f, width = 10, height = 10, units = "in", res = 300)
  tryCatch(print(wrap_plots(pl, ncol = 1)), finally = dev.off())
  pngs = c(pngs, f)
}
stopifnot(length(n_sp) == 26, n_sp[["MR1"]] == 3799)
image_write(image_read(pngs), "final/FigureS4.pdf", format = "pdf",
            density = "300x300", compression = "JPEG", quality = 92)

# ---- 6. Spatial: lesion-level analyses (25 lesions) --------------------------
# 6a. Lesion-level sensitivity analysis of the within-section correlations
sp = fread("results/spatial_within_sample_cor.csv")
sp[, sample  := sub("^GSM[0-9]+_", "", prefix)]
sp[, les_old := sub("[AB]$", "", sample)]                 # original definition (24; merged MR10A/B in error)
sp[, les_new := sub("^MR11[AB]$", "MR11", sample)]        # corrected definition (25)
stopifnot(uniqueN(sp$les_old) == 24, uniqueN(sp$les_new) == 25)
cnt_les = function(col) sp[, .(rho = mean(rho)), by = c(col, "cell", "mn")][,
            .(n = .N, n_pos = sum(rho > 0), median = round(median(rho), 3)), by = .(cell, mn)]
print(merge(cnt_les("les_old"), cnt_les("les_new"), by = c("cell", "mn"), suffixes = c("_old24", "_new25")))

# 6b. Tumour-rich spots per section (tumor_one(), main script lines 1864-1872)
stopifnot(grepl("^tumor_one *= *function", main[1864]))
eval(parse(text = paste(main[1864:1872], collapse = "\n")))
meta = data.table(prefix = sp_pre, title = sec)
meta[, grp2    := fifelse(grepl("^PT", title), "Primary", "Relapse_Met")]
meta[, les_old := sub("[AB]$", "", title)]
meta[, les_new := sub("^MR11[AB]$", "MR11", title)]
sp_tum = merge(rbindlist(lapply(meta$prefix, tumor_one)), meta, by = "prefix")
fwrite(sp_tum, "results/spatial_tumor_spots_by_section.csv")

# 6c. Linear mixed models (patient as random effect), as 19_revision_checks.R lines 1321-1327
tiers = names(mn_sets)
run_lme = function(les) {
  d = sp_tum[, lapply(.SD, mean), by = c(les, "grp2"), .SDcols = tiers]
  setnames(d, les, "lesion")
  d[, pid := paste0("P", gsub("[^0-9]", "", lesion))]
  res = rbindlist(lapply(tiers, function(t) {
    tt = summary(lme(as.formula(paste(t, "~ grp2")), random = ~ 1 | pid, data = d, method = "REML"))$tTable
    data.table(tier = t, diff = round(tt[2, "Value"], 3), p = signif(tt[2, "p-value"], 3))
  }))
  res[, p_adj := signif(p.adjust(p, method = "BH"), 3)]
  list(data = d, res = res)
}
old = run_lme("les_old")                                   # reproduces results/spatial_primary_vs_relapse_lme.csv
new = run_lme("les_new")
print(old$res); print(new$res)
fwrite(new$data, "results/spatial_tumor_spots_by_lesion25.csv")
fwrite(new$res,  "results/spatial_primary_vs_relapse_lme_25lesions.csv")

# 6d. Sensitivity: exclude the patient with two separate relapse lesions (MR10A, MR10B)
fit10 = lme(T3_Glyco ~ grp2, random = ~ 1 | pid, data = new$data[pid != "P10"], method = "REML")
print(summary(fit10)$tTable)

# ---- 7. Fig. 5: final figure (panels A-C unchanged; D with 25 lesions) ------
f5 = readRDS("results/figures/final/fig5_data_19pt.rds")
tier_lab = c(T1_Transport = "T1 transport", T2_Enzyme = "T2 Mn enzymes",
             T3_Glyco = "T3 glycosylation", T4_STING = "T4 cGAS\u2013STING")
les = new$data
td  = melt(les, id.vars = c("lesion", "grp2", "pid"), measure.vars = names(tier_lab),
           variable.name = "tier", value.name = "score")
td[, tier := factor(tier_lab[as.character(tier)], levels = tier_lab)]
td[, grp2 := factor(grp2, levels = c("Primary", "Relapse_Met"))]
pw = data.table(tier = factor(tier_lab[new$res$tier], levels = tier_lab), p = new$res$p)
pw = merge(pw, td[, .(y = max(score)), by = tier], by = "tier")
pw[, lab := sprintf("P = %.3f", p)]
paired = td[, uniqueN(grp2), by = pid][V1 == 2, pid]

# Panel A (19_revision_checks.R lines 1495-1545 and 1575-1581)
a_vars = c(Tumor = "Tumor (osteoblastic)", Myeloid = "Myeloid", T2 = "Mn enzymes",
           Osteoclast = "Osteoclast", aVb3 = "\u03b1v\u03b23 integrin")
one_map = function(v) {
  d = f5$sc5[, .(x, y, val = get(v))]
  ggplot(d, aes(y, x, colour = val)) +
    rasterise(geom_point(size = 0.6, stroke = 0), dpi = 300) +
    scale_colour_distiller(palette = "Spectral", limits = quantile(d$val, c(0.01, 0.99)), oob = scales::squish,
                           name = NULL, breaks = scales::breaks_extended(3),
                           labels = scales::label_number(accuracy = 0.1)) +
    scale_y_reverse() + coord_fixed() + labs(title = a_vars[[v]]) +
    theme_void(base_size = 8, base_family = "Arial") +
    theme(plot.title = element_text(face = "bold", size = 8, hjust = 0.5), legend.position = "bottom",
          legend.key.width = unit(6, "mm"), legend.key.height = unit(1.5, "mm"),
          legend.text = element_text(size = 8), plot.margin = margin(2, 4, 2, 4, "pt"))
}
p5a = wrap_plots(lapply(names(a_vars), one_map), nrow = 1)

# Panel B (lines 1412-1428 and 1473-1478)
hm = copy(f5$hm5)
hm[, cell := factor(cell, levels = rev(c("Myeloid", "Endothelial", "Osteoclast", "Fibroblast/MSC", "T/NK", "Osteoblastic OS")))]
hm[, txt  := sub("^-0\\.00$", "0.00", sprintf("%.2f", rho))]
p5b = ggplot(hm, aes(mn, cell, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.5) + geom_text(aes(label = txt), size = 8 / .pt) +
  scale_fill_gradient2(low = "#4575B4", mid = "white", high = "#B2182B", midpoint = 0,
                       limits = c(-0.15, 0.15), breaks = c(-0.1, 0, 0.1), oob = scales::squish,
                       name = "Median\nSpearman \u03c1") +
  labs(x = NULL, y = NULL) + theme_jtm() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.line = element_blank(),
        axis.ticks = element_blank(), legend.position = "right") +
  guides(fill = guide_colourbar(barwidth = unit(2.5, "mm"), barheight = unit(22, "mm")))

# Panel C (lines 1431-1449, 1481-1488 and 1596)
pc = copy(f5$pc); pos = copy(f5$pos)
lv = sapply(c("Myeloid", "Osteoclast.*\u03b1v\u03b23", "Endothelial", "Fibroblast", "Osteoclast.*transport", "T/NK"),
            function(o) grep(o, levels(factor(f5$pc$lab)), value = TRUE))
pc[,  lab := factor(as.character(lab), levels = rev(lv))]
pos[, lab := factor(as.character(lab), levels = rev(lv))]
set.seed(1)
p5c = ggplot(pc, aes(rho, lab)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  geom_jitter(height = 0.15, width = 0, size = 1.2, colour = "#3C6EA8", alpha = 0.75) +
  geom_errorbar(data = pos, aes(y = lab, xmin = med, xmax = med), inherit.aes = FALSE,
                width = 0.45, colour = "#C0392B", linewidth = 0.6) +
  geom_text(data = pos, aes(x = Inf, y = lab, label = paste0(k, "/", n)), inherit.aes = FALSE,
            hjust = -0.3, size = 8 / .pt) +
  annotate("text", x = Inf, y = Inf, label = "\u03c1 > 0", hjust = -0.3, vjust = -0.6, fontface = "bold", size = 8 / .pt) +
  scale_y_discrete(labels = function(x) gsub("\\s[\u2013\u2014-]\\s", "\n", x)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5), labels = c("0", "0.25", "0.5")) +
  coord_cartesian(clip = "off") + labs(x = "Spearman \u03c1 per patient", y = NULL) + theme_jtm() +
  theme(plot.margin = margin(12, 30, 5, 5, "pt"))

# Panel D (lines 1457-1467, data replaced by 25 lesions)
p5d = ggplot(td, aes(grp2, score, colour = grp2)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, linewidth = 0.5, fill = NA) +
  geom_line(data = td[pid %in% paired], aes(group = pid), colour = "grey70", linewidth = 0.3) +
  geom_point(size = 1.2) +
  geom_text(data = pw, aes(x = 1.5, y = y * 1.08, label = lab), inherit.aes = FALSE, size = 8 / .pt) +
  facet_wrap(~ tier, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = c("#0072B2", "#CC79A7"), guide = "none") +
  scale_x_discrete(labels = c("Primary", "Relapse/\nmetastasis")) +
  labs(x = NULL, y = "Mean score in\ntumor-rich spots") + theme_jtm() +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold", size = 8))

p5bc = (p5b | wrap_elements(full = p5c)) + plot_layout(widths = c(1.15, 1))
fig5 = wrap_elements(full = p5a) / p5bc / p5d +
  plot_layout(heights = c(1, 1.25, 1.1)) + plot_annotation(tag_levels = "A")
save_jtm(fig5, "Figure5", h = 230)

# ---- 8. Additional file 1: Figures S1-S5 in one PDF --------------------------
fs = c("final/FigureS1.pdf", "final/FigureS2.pdf", "final/FigureS3.pdf", "final/FigureS4.pdf",
       "results/figures/final/FigS5_survival_forest.pdf")
stopifnot(all(file.exists(fs)))
pdf_combine(fs, output = "final/Additional_file_1.pdf")
stopifnot(pdf_length("final/Additional_file_1.pdf") == 11)
