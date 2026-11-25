# [local-only, disabled for public release] file.rename("scripts/15_figure_theme.R",
# [local-only, disabled for public release]             "scripts/all_lessons_main_copy.R")
# 416. 列出项目里的所有 R 脚本 ==============================
list.files(".", pattern = "\\.R$", recursive = TRUE, full.names = TRUE)
length(readLines("scripts/15_figure_theme.R"))     # 应约 75 行
# 417. 主脚本改名并备份 ==============================
# [local-only, disabled for public release] file.rename("scripts/all_lessons_main_copy.R",
# [local-only, disabled for public release]             "scripts/00_main_all_lessons.R")
file.copy("scripts/00_main_all_lessons.R",
          paste0("scripts/00_main_backup_", format(Sys.Date(), "%Y%m%d"), ".R"))
# 418. 用代码生成主题脚本 ==============================
writeLines(r"---(
# ===== 15_figure_theme.R：JTM 统一作图规格 =====
library(ggplot2)
library(patchwork)
library(data.table)
library(ggrepel)

W_FULL    = 176
W_HALF    = 85
final_dir = "results/figures/final"
dir.create(final_dir, showWarnings = FALSE, recursive = TRUE)

save_jtm = function(p, name, h, w = W_FULL) {
  ggsave(file.path(final_dir, paste0(name, ".tiff")), p, width = w, height = h,
         units = "mm", dpi = 300, device = ragg::agg_tiff, compression = "lzw")
  ggsave(file.path(final_dir, paste0(name, ".pdf")),  p, width = w, height = h,
         units = "mm", device = cairo_pdf)
}

theme_jtm = function() {
  theme_classic(base_size = 8, base_family = "Arial") +
  theme(axis.text        = element_text(size = 8, colour = "black"),
        axis.title       = element_text(size = 9),
        legend.text      = element_text(size = 8),
        legend.title     = element_text(size = 8),
        strip.text       = element_text(size = 8, face = "bold"),
        strip.background = element_blank(),
        plot.title       = element_text(size = 9, face = "bold"),
        axis.line        = element_line(linewidth = 0.3),
        axis.ticks       = element_line(linewidth = 0.3),
        plot.tag         = element_text(size = 12, face = "bold"))
}
theme_set(theme_jtm())

ct_lab = c(OS_Osteoblastic   = "Osteoblastic OS",  OS_Chondroblastic = "Chondroblastic OS",
           Fibroblast_MSC    = "Fibroblast/MSC",   Myeloid           = "Myeloid",
           Osteoclast        = "Osteoclast",       Endothelial       = "Endothelial",
           Pericyte          = "Pericyte",         T_NK              = "T/NK",
           B_Plasma          = "B/plasma",         Mast              = "Mast",
           Lung_epithelial   = "Lung epithelial",  Myocyte           = "Myocyte")
tier_lab = c(T1_Transport = "T1 transport",     T2_Enzyme = "T2 Mn enzymes",
             T3_Glyco     = "T3 glycosylation", T4_STING  = "T4 cGAS–STING")

ct_col = c("Osteoblastic OS"   = "#3B75AF", "Chondroblastic OS" = "#8FB8DE",
           "Fibroblast/MSC"    = "#E69F00", "Myeloid"           = "#009E73",
           "Osteoclast"        = "#CC79A7", "Endothelial"       = "#D55E00",
           "Pericyte"          = "#B8860B", "T/NK"              = "#7B3294",
           "B/plasma"          = "#C2A5CF", "Mast"              = "#999999",
           "Lung epithelial"   = "#66C2A5", "Myocyte"           = "#4D4D4D")

fmt_p = function(p) ifelse(p < 0.001, sprintf("%.1e", p), sprintf("%.3f", p))

plot_umap = function(obj) {
  d   = as.data.table(Seurat::Embeddings(obj, "umap"))
  setnames(d, c("UMAP1", "UMAP2"))
  d$ct = obj$ct
  lab = d[, .(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2)), by = ct]
  ggplot(d, aes(UMAP1, UMAP2, colour = ct)) +
    ggrastr::rasterise(geom_point(size = 0.05, stroke = 0), dpi = 300) +
    geom_text_repel(data = lab, aes(label = ct), colour = "black",
                    bg.colour = "white", bg.r = 0.15,
                    size = 8 / .pt, min.segment.length = Inf, seed = 1) +
    scale_colour_manual(values = ct_col) +
    labs(x = "UMAP 1", y = "UMAP 2") +
    theme(legend.position = "none")
}
)---", "scripts/15_figure_theme.R", useBytes = TRUE)

length(readLines("scripts/15_figure_theme.R", encoding = "UTF-8"))
grep("install", readLines("scripts/15_figure_theme.R"), value = TRUE)
source("scripts/15_figure_theme.R", encoding = "UTF-8")
library(CellChat)
# 408. 读入 CellChat 对象并提取通讯表 ==============================
cellchat = readRDS("data_proc/cellchat.rds")
cc_df    = as.data.table(subsetCommunication(cellchat))
print(dim(cc_df))
print(colnames(cc_df))
print(levels(cellchat@idents))

# 409. 查看关键配体–受体对 ==============================
key = cc_df[grepl("MIF|^COL|SPP1|ITGA10|CD74", interaction_name)]
print(key[order(-prob)][1:40, .(source, target, interaction_name_2, prob, pval)])

# 410. 查看细胞间通讯数量与强度矩阵 ==============================
print(round(cellchat@net$count, 0))
print(round(cellchat@net$weight, 2))
# 419. 定义细胞顺序与矩阵转换函数 ==============================
ct_order = c("Myeloid", "Endothelial", "Pericyte", "Osteoclast", "Fibroblast/MSC",
             "T/NK", "Chondroblastic OS", "Osteoblastic OS")
mat_long = function(m) {
  d = as.data.table(as.table(m))
  setnames(d, c("source", "target", "value"))
  d[, source := factor(unname(ct_lab[source]), levels = rev(ct_order))]
  d[, target := factor(unname(ct_lab[target]), levels = ct_order)]
  d
}
plot_mat = function(d, title) {
  ggplot(d, aes(target, source, fill = value)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    scale_fill_gradient(low = "white", high = "#B2182B", name = "Strength") +
    labs(x = "Receiver", y = "Sender", title = title) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.line   = element_blank(),
          axis.ticks  = element_blank())
}

# 420. 重画 Fig 4A：全部通讯强度 ==============================
p4a = plot_mat(mat_long(cellchat@net$weight), "All interactions")

# 421. 重画 Fig 4B：整合素通讯强度 ==============================
itg = cc_df[grepl("ITG", receptor), .(value = sum(prob)), by = .(source, target)]
itg = merge(CJ(source = levels(cellchat@idents), target = levels(cellchat@idents)),
            itg, by = c("source", "target"), all.x = TRUE)
itg[is.na(value), value := 0]
itg[, source := factor(unname(ct_lab[source]), levels = rev(ct_order))]
itg[, target := factor(unname(ct_lab[target]), levels = ct_order)]
p4b = plot_mat(itg, "Integrin receptors") +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank())

# 422. 重画 Fig 4C：精选配体–受体气泡图 ==============================
pairs = c("MIF - (CD74+CXCR4)",      "MIF - (CD74+CD44)",       "APP - CD74",
          "COL1A1 - (ITGA10+ITGB1)", "COL1A2 - (ITGA10+ITGB1)", "COL1A1 - (ITGA1+ITGB1)",
          "COL1A1 - (ITGA11+ITGB1)", "SPP1 - (ITGAV+ITGB1)",    "SPP1 - (ITGA5+ITGB1)",
          "SPP1 - CD44")
sel = cc_df[trimws(interaction_name_2) %in% pairs]
sel[, pair := factor(trimws(interaction_name_2), levels = rev(pairs))]
sel[, link := paste0(ct_lab[as.character(source)], " → ", ct_lab[as.character(target)])]
must      = c("Chondroblastic OS → Myeloid", "Endothelial → Myeloid",
              "Fibroblast/MSC → Chondroblastic OS", "Osteoclast → Osteoclast")
ranked    = sel[, .(s = sum(prob)), by = link][order(-s), link]
top_links = unique(c(must, ranked))[1:14]
print(top_links)
sel = sel[link %in% top_links]
sel[, link := factor(link, levels = top_links)]
p4c = ggplot(sel, aes(link, pair)) +
  geom_point(aes(size = prob, colour = prob)) +
  scale_colour_gradient(low = "#9ECAE1", high = "#08306B", name = "Communication\nprobability") +
  scale_size_continuous(range = c(1, 4), guide = "none") +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x      = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_line(linewidth = 0.2, colour = "grey92"))

# 423. 重画 Fig 4D：受体表达点图 ==============================
sc    = readRDS("data_proc/sc_tme_scored.rds")
sc$ct = factor(unname(ct_lab[as.character(sc$celltype)]), levels = rev(ct_order))
rec   = c("CD74", "CXCR4", "CD44", "ITGB1", "ITGA1", "ITGA10", "ITGA11", "ITGA5", "ITGAV", "ITGB3")
p4d = Seurat::DotPlot(sc, features = rec, group.by = "ct",
                      cols = c("grey85", "#3B75AF"), dot.scale = 3.5) +
  theme_jtm() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic")) +
  guides(colour = guide_colourbar(title = "Scaled mean\nexpression"),
         size   = guide_legend(title = "% expressing")) +
  labs(x = NULL, y = NULL)

# 424. 拼合并导出
file.info(file.path(final_dir, "Figure4.tiff"))[, c("size", "mtime")]
sapply(c("p4a", "p4b", "p4c", "sc", "p4d", "fig4"), exists)
# 424. 拼合并导出 Fig 4 ==============================
fig4 = (p4a | p4b) / p4c / p4d +
  plot_layout(heights = c(1, 1.1, 0.8)) +
  plot_annotation(tag_levels = "A")
save_jtm(fig4, "Figure4", h = 220)
saveRDS(list(weight = cellchat@net$weight, itg = itg, sel = sel),
        "results/figures/final/fig4_data.rds")
file.info(file.path(final_dir, "Figure4.tiff"))[, c("size", "mtime")]
# 425. 修正 Fig 4D 图例为横向并排 ==============================
p4d = p4d +
  guides(colour = guide_colourbar(title = "Scaled mean expression", title.position = "top",
                                  direction = "horizontal",
                                  barwidth = unit(22, "mm"), barheight = unit(2.5, "mm")),
         size   = guide_legend(title = "% expressing", title.position = "top", nrow = 1)) +
  theme(legend.position = "right",
        legend.box      = "vertical",
        legend.margin   = margin(0, 0, 0, 0))

fig4 = (p4a | p4b) / p4c / p4d +
  plot_layout(heights = c(1, 1.1, 0.85)) +
  plot_annotation(tag_levels = "A")
save_jtm(fig4, "Figure4", h = 225)
# 426. 查看 Fig 5 所需数据 ==============================
sp_cor = fread("results/spatial_within_sample_cor.csv")
sp_all = fread("results/spatial_all_pairs_cor.csv")
sp_tum = fread("results/spatial_tumor_spots_by_patient.csv")
for (d in list(sp_cor, sp_all, sp_tum)) { print(colnames(d)); print(head(d, 3)) }
print(unique(sp_cor[[2]]))

# 427. 查看空间样本文件夹 ==============================
list.dirs("data_raw/spatial/samples", recursive = FALSE)
list.files("data_proc", pattern = "spatial|visium", ignore.case = TRUE)
# 428. 空间数据的标签对照 ==============================
sp_cell_lab = c(ct_lab, Fibroblast = "Fibroblast/MSC", Tumour = "Tumor", Tumor = "Tumor")
mn_lab      = c(T1_Transport = "Mn transport", T2_Enzyme = "Mn enzymes",
                T3_Glyco     = "Mn glycosylation", T4_STING = "cGAS–STING", aVb3 = "αvβ3 integrin")
fix_lab     = function(x, map) ifelse(x %in% names(map), map[x], gsub("_", " ", x))
print(unique(sp_all$cell)); print(unique(sp_all$mn))

# 429. 重画 Fig 5B：组织内相关系数中位数热图 ==============================
hm5 = sp_all[, .(rho = median(rho)), by = .(cell, mn)]
hm5[, cell := fix_lab(cell, sp_cell_lab)]
hm5[, mn   := factor(fix_lab(mn, mn_lab), levels = mn_lab)]
p5b = ggplot(hm5, aes(mn, cell, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 8 / .pt) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       name = "Median\nSpearman ρ") +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line   = element_blank(),
        axis.ticks  = element_blank())

# 430. 重画 Fig 5C：每位病人的组织内相关系数 ==============================
pc = sp_cor[, .(rho = mean(rho)), by = .(patient, cell, mn)]
pc[, lab := paste0(fix_lab(cell, sp_cell_lab), " — ", fix_lab(mn, mn_lab))]
pos = pc[, .(k = sum(rho > 0), n = .N, med = median(rho)), by = lab][order(-med)]
pc[, lab := factor(lab, levels = rev(pos$lab))]
pos[, lab := factor(lab, levels = rev(pos$lab))]
print(pos)
p5c = ggplot(pc, aes(rho, lab)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3, colour = "grey40") +
  geom_point(size = 0.9, colour = "#3B75AF", alpha = 0.7,
             position = position_jitter(height = 0.15, seed = 1)) +
  geom_point(data = pos, aes(med, lab), shape = 124, size = 4, colour = "#B2182B") +
  geom_text(data = pos, aes(x = Inf, y = lab, label = paste0(k, "/", n)),
            hjust = -0.2, size = 8 / .pt) +
  annotate("text", x = Inf, y = Inf, label = "ρ > 0", hjust = -0.2, vjust = -0.6,
           size = 8 / .pt, fontface = "bold") +
  coord_cartesian(clip = "off") +
  labs(x = "Within-section Spearman ρ (per patient)", y = NULL) +
  theme(plot.margin = margin(15, 30, 5, 5))

# 431. 重画 Fig 5D：原发与复发/转移的肿瘤区评分 ==============================
td = melt(sp_tum, id.vars = c("patient", "group", "grp2"),
          measure.vars = c("T1_Transport", "T2_Enzyme", "T3_Glyco", "T4_STING"),
          variable.name = "tier", value.name = "score")
td[, tier := factor(tier_lab[as.character(tier)], levels = tier_lab)]
td[, grp2 := factor(ifelse(grp2 == "Primary", "Primary", "Relapse/metastasis"),
                    levels = c("Primary", "Relapse/metastasis"))]
pw = td[, .(p = wilcox.test(score ~ grp2)$p.value, y = max(score)), by = tier]
print(pw); print(table(sp_tum$grp2))
p5d = ggplot(td, aes(grp2, score, colour = grp2)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.3, width = 0.6) +
  geom_point(size = 0.9, position = position_jitter(width = 0.12, seed = 1)) +
  geom_text(data = pw, inherit.aes = FALSE, aes(x = 1.5, y = y * 1.06, label = paste0("P = ", fmt_p(p))),
            size = 8 / .pt) +
  facet_wrap(~ tier, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = c("Primary" = "#0072B2", "Relapse/metastasis" = "#CC79A7")) +
  labs(x = NULL, y = "Mean score in\ntumor-rich spots") +
  theme(legend.position = "none")

# 432. 在主脚本中查找 Fig 5A 的代码位置 ==============================
main = readLines("scripts/00_main_all_lessons.R", encoding = "UTF-8")
hits = grep("Load10X_Spatial|Tumour_s|aVb3_s|smooth", main)
print(hits)
print(main[hits])
# 433. 打印 Fig 5A 原始代码 ==============================
cat(main[1895:1950], sep = "\n")
cat(main[1895:1922], sep = "\n")
cat(main[1923:1950], sep = "\n")
# 434. 查找基因集和路径的定义 ==============================
idx = grep("^\\s*(cell_sets|mn_sets|sp_dir)\\s*=", main)
print(idx)
for (i in idx) { cat("---- line", i, "----\n"); cat(main[i:(i + 12)], sep = "\n") }
# 435. 定义空间基因集并读入代表性样本 ==============================
library(Seurat)
sp_dir    = "data_raw/spatial"
cell_sets = list(OS_Osteoblastic = c("RUNX2", "ALPL", "IBSP", "COL1A1"),
                 Osteoclast      = c("ACP5", "CTSK", "MMP9", "DCSTAMP"),
                 Myeloid         = c("CD74", "LYZ", "CD68", "CD14"))
t2_genes  = c("SOD2", "ARG2", "GLUL", "PC", "PPM1A", "PPM1B")
rep_cand  = sp_all[cell == "Myeloid" & mn == "T2_Enzyme" & n_spots > 1000]
rep_cand[, dist := abs(rho - median(rho))]
rep_p     = rep_cand[order(dist)][1, prefix]
print(rep_cand[order(dist)][1:3, .(prefix, rho, n_spots)])
st        = NormalizeData(Load10X_Spatial(file.path(sp_dir, "samples", rep_p)))

# 436. 计算五个展示评分并做邻域平滑 ==============================
show_sets = list(Tumor = cell_sets$OS_Osteoblastic, Myeloid = cell_sets$Myeloid,
                 T2    = t2_genes,                   Osteoclast = cell_sets$Osteoclast)
show_sets = lapply(show_sets, intersect, rownames(st))
st        = AddModuleScore(st, features = show_sets, name = "Show")
ex        = GetAssayData(st, layer = "data")
sc5       = data.table(Tumor      = st$Show1, Myeloid = st$Show2,
                       T2         = st$Show3, Osteoclast = st$Show4,
                       aVb3       = colMeans(as.matrix(ex[intersect(c("ITGAV", "ITGB3"), rownames(ex)), , drop = FALSE])))
coords    = GetTissueCoordinates(st)[colnames(st), ]
D         = as.matrix(dist(coords[, c("x", "y")]))
nn        = t(apply(D, 1, function(r) order(r)[1:7]))
rm(D); gc()
smooth    = function(v) rowMeans(matrix(v[nn], nrow = nrow(nn)))
sc5       = sc5[, lapply(.SD, smooth)]
sc5[, `:=`(x = coords$x, y = coords$y)]

# 437. 重画 Fig 5A：五个并排的空间分布图 ==============================
feat_title = c(Tumor = "Tumor (osteoblastic)", Myeloid = "Myeloid",
               T2 = "Mn enzymes", Osteoclast = "Osteoclast", aVb3 = "αvβ3 integrin")
sp_panel = function(f) {
  ggplot(sc5, aes(y, x, colour = .data[[f]])) +
    ggrastr::rasterise(geom_point(size = 0.35, stroke = 0), dpi = 300) +
    scale_colour_distiller(palette = "Spectral", direction = -1, name = NULL,
                           guide = guide_colourbar(direction = "horizontal",
                                                   barwidth = unit(20, "mm"),
                                                   barheight = unit(2, "mm"))) +
    scale_y_reverse() +
    coord_fixed() +
    labs(title = feat_title[f]) +
    theme_void(base_size = 8, base_family = "Arial") +
    theme(plot.title      = element_text(size = 8, face = "bold", hjust = 0.5),
          legend.position = "bottom",
          legend.text     = element_text(size = 8))
}
p5a = wrap_elements(full = wrap_plots(lapply(names(feat_title), sp_panel), nrow = 1))

# 438. 拼合并导出 Fig 5 ==============================
fig5 = p5a / (p5b | p5c) / p5d +
  plot_layout(heights = c(0.8, 1, 0.75)) +
  plot_annotation(tag_levels = "A")
save_jtm(fig5, "Figure5", h = 215)
saveRDS(list(sc5 = sc5, rep_p = rep_p, hm5 = hm5, pc = pc, pos = pos, td = td, pw = pw),
        "results/figures/final/fig5_data.rds")
# 439. 修正 Fig 5A：色条刻度与顶部留白 ==============================
sp_panel = function(f) {
  ggplot(sc5, aes(y, x, colour = .data[[f]])) +
    ggrastr::rasterise(geom_point(size = 0.35, stroke = 0), dpi = 300) +
    scale_colour_distiller(palette = "Spectral", direction = -1, name = NULL,
                           breaks = scales::breaks_pretty(n = 3),
                           guide  = guide_colourbar(direction = "horizontal",
                                                    barwidth = unit(22, "mm"),
                                                    barheight = unit(2, "mm"))) +
    scale_y_reverse() +
    coord_fixed() +
    labs(title = feat_title[f]) +
    theme_void(base_size = 8, base_family = "Arial") +
    theme(plot.title      = element_text(size = 8, face = "bold", hjust = 0.5),
          legend.position = "bottom",
          legend.text     = element_text(size = 8),
          plot.margin     = margin(12, 2, 2, 2))
}
p5a = wrap_elements(full = wrap_plots(lapply(names(feat_title), sp_panel), nrow = 1))

# 440. 修正 Fig 5B：统一行顺序，图例移到底部 ==============================
sp_order = c("Myeloid", "Endothelial", "Osteoclast", "Fibroblast/MSC", "T/NK", "Osteoblastic OS")
hm5[, cell := factor(cell, levels = rev(sp_order))]
p5b = p5b +
  guides(fill = guide_colourbar(title = "Median Spearman ρ", title.position = "top",
                                direction = "horizontal",
                                barwidth = unit(30, "mm"), barheight = unit(2.5, "mm"))) +
  theme(legend.position = "bottom")

# 441. 修正 Fig 5C 横轴标题与 Fig 5D 分组标签 ==============================
p5c = p5c + labs(x = "Within-section Spearman ρ")
td[, grp2 := factor(ifelse(grp2 == "Primary", "Primary", "Relapse/\nmetastasis"),
                    levels = c("Primary", "Relapse/\nmetastasis"))]
p5d = p5d %+% td +
  scale_colour_manual(values = c("Primary" = "#0072B2", "Relapse/\nmetastasis" = "#CC79A7"))

# 442. 重新拼合并导出 Fig 5 ==============================
fig5 = p5a / (p5b | free(p5c)) / p5d +
  plot_layout(heights = c(0.8, 1, 0.75)) +
  plot_annotation(tag_levels = "A")
fig5[[2]] = fig5[[2]] + plot_layout(widths = c(1.15, 1))
save_jtm(fig5, "Figure5", h = 215)
# 442. 重新拼合并导出 Fig 5（修正版） ==============================
row2 = (p5b | p5c) + plot_layout(widths = c(1.15, 1))
fig5 = p5a / row2 / p5d +
  plot_layout(heights = c(0.8, 1, 0.75)) +
  plot_annotation(tag_levels = "A")
save_jtm(fig5, "Figure5", h = 215)
# 443. 修正 Fig 5B 图例刻度与负零 ==============================
p5b = ggplot(hm5, aes(mn, cell, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", round(rho, 2) + 0)), size = 8 / .pt) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       breaks = c(-0.1, 0, 0.1),
                       guide  = guide_colourbar(title = "Median Spearman ρ", title.position = "top",
                                                direction = "horizontal",
                                                barwidth = unit(30, "mm"), barheight = unit(2.5, "mm"))) +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x     = element_text(angle = 45, hjust = 1),
        axis.line       = element_blank(),
        axis.ticks      = element_blank(),
        legend.position = "bottom")

# 444. C 单独打包，避免横轴标题被拉远 ==============================
row2 = (p5b | wrap_elements(full = p5c)) + plot_layout(widths = c(1.15, 1))
fig5 = p5a / row2 / p5d +
  plot_layout(heights = c(0.8, 1, 0.75)) +
  plot_annotation(tag_levels = "A")
save_jtm(fig5, "Figure5", h = 215)
# 445. 修正 Fig 5C：标签分两行，加宽绘图区 ==============================
pc[,  lab2 := sub(" — ", "\n", as.character(lab))]
pos[, lab2 := sub(" — ", "\n", as.character(lab))]
lv = sub(" — ", "\n", levels(pc$lab))
pc[,  lab2 := factor(lab2, levels = lv)]
pos[, lab2 := factor(lab2, levels = lv)]
p5c = ggplot(pc, aes(rho, lab2)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3, colour = "grey40") +
  geom_point(size = 0.9, colour = "#3B75AF", alpha = 0.7,
             position = position_jitter(height = 0.15, seed = 1)) +
  geom_point(data = pos, aes(med, lab2), shape = 124, size = 4, colour = "#B2182B") +
  geom_text(data = pos, aes(x = Inf, y = lab2, label = paste0(k, "/", n)),
            hjust = -0.2, size = 8 / .pt) +
  annotate("text", x = Inf, y = Inf, label = "ρ > 0", hjust = -0.2, vjust = -0.6,
           size = 8 / .pt, fontface = "bold") +
  scale_x_continuous(breaks = c(0, 0.25, 0.5)) +
  coord_cartesian(clip = "off") +
  labs(x = "Within-section Spearman ρ", y = NULL) +
  theme(plot.margin = margin(15, 30, 5, 5))

# 446. 重新拼合并导出 Fig 5 ==============================
row2 = (p5b | wrap_elements(full = p5c)) + plot_layout(widths = c(1, 1))
fig5 = p5a / row2 / p5d +
  plot_layout(heights = c(0.8, 1, 0.75)) +
  plot_annotation(tag_levels = "A")
save_jtm(fig5, "Figure5", h = 215)
# 447. 用原始列名直接生成两行标签 ==============================
pc[, lab2 := paste0(fix_lab(cell, sp_cell_lab), "\n", fix_lab(mn, mn_lab))]
pos2 = pc[, .(k = sum(rho > 0), n = .N, med = median(rho)), by = lab2][order(-med)]
pc[,   lab2 := factor(lab2, levels = rev(pos2$lab2))]
pos2[, lab2 := factor(lab2, levels = rev(pos2$lab2))]
print(pos2)

p5c = ggplot(pc, aes(rho, lab2)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3, colour = "grey40") +
  geom_point(size = 0.9, colour = "#3B75AF", alpha = 0.7,
             position = position_jitter(height = 0.15, seed = 1)) +
  geom_point(data = pos2, aes(med, lab2), shape = 124, size = 4, colour = "#B2182B") +
  geom_text(data = pos2, aes(x = Inf, y = lab2, label = paste0(k, "/", n)),
            hjust = -0.2, size = 8 / .pt) +
  annotate("text", x = Inf, y = Inf, label = "ρ > 0", hjust = -0.2, vjust = -0.6,
           size = 8 / .pt, fontface = "bold") +
  scale_x_continuous(breaks = c(0, 0.25, 0.5)) +
  coord_cartesian(clip = "off") +
  labs(x = "Within-section Spearman ρ", y = NULL) +
  theme(plot.margin = margin(15, 30, 5, 5))

row2 = (p5b | wrap_elements(full = p5c)) + plot_layout(widths = c(1, 1))
fig5 = p5a / row2 / p5d +
  plot_layout(heights = c(0.8, 1, 0.75)) +
  plot_annotation(tag_levels = "A")
save_jtm(fig5, "Figure5", h = 215)
# 448. 保存 Fig 5 绘图数据 ==============================
saveRDS(list(sc5 = sc5, rep_p = rep_p, hm5 = hm5, pc = pc, pos = pos2, td = td, pw = pw),
        "results/figures/final/fig5_data.rds")

# 449. 查看 Fig 6 所需数据 ==============================
imm  = fread("results/gse21257_immune_scores.csv")
tclin = fread("data_proc/target_clin_scores.csv")
gclin = fread("data_proc/gse21257_clin_scores.csv")
null = readRDS("results/random_null_gse21257.rds")
rnd  = fread("results/random_geneset_gse21257.csv")
for (d in list(imm, tclin, gclin, rnd)) { print(colnames(d)); print(head(d, 2)) }
str(null, max.level = 1)
file.info(file.path("results/figures/final",
                    c("Figure5.tiff", "fig5_data.rds")))[, c("size", "mtime")]
source("scripts/15_figure_theme.R", encoding = "UTF-8")
# 450. 重新计算 TARGET 的免疫细胞评分 ==============================
texpr = readRDS("data_proc/target_expr.rds")
print(dim(texpr)); print(class(texpr))
imm_sets = list(T_cell      = c("CD3D", "CD3E", "CD3G", "CD2", "CD5", "LCK"),
                CD8_T       = c("CD8A", "CD8B"),
                Macrophage  = c("CD68", "CD163", "CD14", "C1QA", "C1QB"),
                Endothelial = c("VWF", "PECAM1", "CDH5", "EMCN"),
                Osteoclast  = c("ACP5", "CTSK", "MMP9"))
tx    = as.matrix(texpr)
z_t   = t(scale(t(tx)))
cs_t  = sapply(imm_sets, function(g) colMeans(z_t[intersect(g, rownames(z_t)), , drop = FALSE]))
sting_t = tclin$T4_STING[match(colnames(tx), tclin$sample)]
print(sum(!is.na(sting_t)))                                # 应为 88

# 451. 两个队列的 STING–细胞相关系数（核对原图） ==============================
cor_t = apply(cs_t, 2, cor, y = sting_t, method = "spearman", use = "complete.obs")
cor_g = sapply(names(imm_sets), function(k) cor(imm$STING, imm[[k]], method = "spearman"))
print(round(rbind(TARGET = cor_t, GSE21257 = cor_g), 2))

# 452. 逻辑回归：转移 ~ STING（± 免疫细胞） ==============================
d6  = imm[, .(meta  = as.integer(meta),
              STING = scale(STING)[, 1],
              Tcell = scale(T_cell)[, 1],
              Mac   = scale(Macrophage)[, 1])]
fit = function(f, model) {
  m  = glm(f, data = d6, family = binomial)
  ci = confint.default(m)
  data.table(model = model, term = names(coef(m))[-1],
             OR = exp(coef(m)[-1]), lo = exp(ci[-1, 1]), hi = exp(ci[-1, 2]),
             p  = summary(m)$coefficients[-1, 4])
}
lr = rbind(fit(meta ~ STING,         "STING alone"),
           fit(meta ~ STING + Tcell, "STING + T cells"),
           fit(meta ~ STING + Mac,   "STING + macrophages"))
print(lr); print(table(imm$meta))
# 453. 重画 Fig 6A：STING 与免疫细胞评分相关性热图 ==============================
imm_lab = c(T_cell = "T cells", CD8_T = "CD8+ T cells", Macrophage = "Macrophages",
            Endothelial = "Endothelial cells", Osteoclast = "Osteoclasts")
coh_lab = c(TARGET = "TARGET-OS\n(RNA-seq, n = 88)", GSE21257 = "GSE21257\n(array, n = 53)")
h6 = rbind(data.table(cohort = "TARGET",   cell = names(cor_t), rho = cor_t),
           data.table(cohort = "GSE21257", cell = names(cor_g), rho = cor_g))
h6[, cohort := factor(coh_lab[cohort], levels = coh_lab)]
h6[, cell   := factor(imm_lab[cell],   levels = rev(imm_lab))]
p6a = ggplot(h6, aes(cohort, cell, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 8 / .pt) +
  scale_fill_gradient(low = "white", high = "#B2182B", limits = c(0, 0.8),
                      guide = guide_colourbar(title = "Spearman ρ with\nSTING score",
                                              barheight = unit(25, "mm"), barwidth = unit(2.5, "mm"))) +
  labs(x = NULL, y = NULL) +
  theme(axis.line = element_blank(), axis.ticks = element_blank())

# 454. 重画 Fig 6B：转移的逻辑回归森林图 ==============================
term_lab = c(STING = "STING score", Tcell = "T-cell score", Mac = "Macrophage score")
mod_lv   = c("STING alone", "STING + T cells", "STING + macrophages")
lr[, y    := match(model, mod_lv) + ifelse(term == "STING", -0.15, 0.15)]
lr[, term := factor(term_lab[term], levels = term_lab)]
p6b = ggplot(lr, aes(OR, y, colour = term)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.3, colour = "grey40") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.4) +
  geom_point(size = 1.6) +
  scale_x_log10(breaks = c(0.1, 0.3, 1, 2)) +
  scale_y_reverse(breaks = 1:3, labels = mod_lv) +
  scale_colour_manual(values = c("STING score" = "black", "T-cell score" = "grey55",
                                 "Macrophage score" = "#009E73"), name = NULL) +
  labs(x = "OR for metastasis per SD (log scale)", y = NULL,
       title = "GSE21257: metastasis (34 vs 19)") +
  theme(legend.position = "bottom")

# 455. 重画 Fig 6C：随机基因集零分布 ==============================
fmt_p2 = function(p) formatC(signif(p, 2), format = "fg")
test_lab = c(STING_meta = "STING score vs metastasis",
             Glyco_cox  = "Mn glycosylation vs overall survival")
nd = rbind(data.table(test = test_lab["STING_meta"], lp = -log10(null$null_meta)),
           data.table(test = test_lab["Glyco_cox"],  lp = -log10(null$null_cox)))
od = data.table(test = test_lab[names(null$obs)], lp = -log10(null$obs),
                lab  = paste0("Observed P = ", fmt_p2(null$obs),
                              "\nEmpirical P = ", fmt_p2(null$emp)))
nd[, test := factor(test, levels = test_lab)]
od[, test := factor(test, levels = test_lab)]
p6c = ggplot(nd, aes(lp)) +
  geom_histogram(bins = 40, fill = "grey75", colour = NA) +
  geom_vline(data = od, aes(xintercept = lp), colour = "#B2182B", linewidth = 0.6) +
  geom_text(data = od, aes(x = lp, y = Inf, label = lab), hjust = -0.05, vjust = 1.3,
            size = 8 / .pt, colour = "#B2182B", lineheight = 0.9) +
  facet_wrap(~ test, nrow = 1, scales = "free") +
  labs(x = expression(-log[10](P) ~ "of 1,000 random gene sets of equal size"), y = "Count")

# 456. 拼合并导出 Fig 6 ==============================
fig6 = (p6a | p6b) / p6c +
  plot_layout(heights = c(1, 0.9)) +
  plot_annotation(tag_levels = "A")
fig6[[1]] = NULL
fig6 = ((p6a | p6b) + plot_layout(widths = c(0.8, 1.2))) / p6c +
  plot_layout(heights = c(1, 0.9)) +
  plot_annotation(tag_levels = "A")
save_jtm(fig6, "Figure6", h = 150)
saveRDS(list(h6 = h6, lr = lr, nd = nd, od = od),
        "results/figures/final/fig6_data.rds")
# 457. 修正 Fig 6A：队列名称分三行 ==============================
coh_lab3 = c("TARGET-OS\n(RNA-seq,\nn = 88)", "GSE21257\n(array,\nn = 53)")
h6[, cohort := factor(coh_lab3[as.integer(cohort)], levels = coh_lab3)]
p6a = p6a %+% h6

# 458. 重画 Fig 6B：标准森林图（每行一个系数） ==============================
row_lab = c("STING alone|STING score"              = "STING (alone)",
            "STING + T cells|STING score"          = "STING (+ T cells)",
            "STING + T cells|T-cell score"         = "T cells",
            "STING + macrophages|STING score"      = "STING (+ macrophages)",
            "STING + macrophages|Macrophage score" = "Macrophages")
lr[, row := factor(row_lab[paste(model, term, sep = "|")], levels = rev(row_lab))]
lr[, txt := sprintf("%.2f (%.2f–%.2f)", OR, lo, hi)]
p6b = ggplot(lr, aes(OR, row, colour = term)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.3, colour = "grey40") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.4) +
  geom_point(size = 1.6) +
  geom_text(aes(x = Inf, label = txt), hjust = -0.1, size = 8 / .pt, colour = "black") +
  annotate("text", x = Inf, y = Inf, label = "OR (95% CI)", hjust = -0.1, vjust = -0.6,
           size = 8 / .pt, fontface = "bold") +
  scale_x_log10(breaks = c(0.1, 0.3, 1, 2)) +
  scale_colour_manual(values = c("STING score" = "black", "T-cell score" = "grey55",
                                 "Macrophage score" = "#009E73")) +
  coord_cartesian(clip = "off") +
  labs(x = "OR for metastasis per SD (log scale)", y = NULL,
       title = "GSE21257: metastasis (34 vs 19)") +
  theme(legend.position = "none", plot.margin = margin(15, 70, 5, 5))

# 459. 修正 Fig 6C：P 值文字放在竖线左侧 ==============================
p6c = ggplot(nd, aes(lp)) +
  geom_histogram(bins = 40, fill = "grey75", colour = NA) +
  geom_vline(data = od, aes(xintercept = lp), colour = "#B2182B", linewidth = 0.6) +
  geom_text(data = od, aes(x = lp, y = Inf, label = lab), hjust = 1.05, vjust = 1.3,
            size = 8 / .pt, colour = "#B2182B", lineheight = 0.9) +
  facet_wrap(~ test, nrow = 1, scales = "free") +
  labs(x = expression(-log[10](P) ~ "of 1,000 random gene sets of equal size"), y = "Count")

# 460. 重新拼合并导出 Fig 6 ==============================
fig6 = ((p6a | p6b) + plot_layout(widths = c(1, 1.3))) / p6c +
  plot_layout(heights = c(1, 0.9)) +
  plot_annotation(tag_levels = "A")
save_jtm(fig6, "Figure6", h = 150)
saveRDS(list(h6 = h6, lr = lr, nd = nd, od = od),
        "results/figures/final/fig6_data.rds")
# 461. 修正 Fig 6A 横轴标签 ==============================
coh_lab4 = c("TARGET-OS\n(n = 88)", "GSE21257\n(n = 53)")
h6[, cohort := factor(coh_lab4[as.integer(cohort)], levels = coh_lab4)]
p6a = p6a %+% h6

# 462. 修正 Fig 6B 标题位置 ==============================
p6b = p6b + theme(plot.title.position = "plot")

# 463. 修正 Fig 6C：左格文字放右侧，右格放左侧 ==============================
od[, hj := ifelse(test == test_lab["STING_meta"], -0.05, 1.05)]
p6c = ggplot(nd, aes(lp)) +
  geom_histogram(bins = 40, fill = "grey75", colour = NA) +
  geom_vline(data = od, aes(xintercept = lp), colour = "#B2182B", linewidth = 0.6) +
  geom_text(data = od, aes(x = lp, y = Inf, label = lab, hjust = hj), vjust = 1.3,
            size = 8 / .pt, colour = "#B2182B", lineheight = 0.9) +
  facet_wrap(~ test, nrow = 1, scales = "free") +
  labs(x = expression(-log[10](P) ~ "of 1,000 random gene sets of equal size"), y = "Count")

# 464. 重新拼合并导出 Fig 6 ==============================
fig6 = ((p6a | p6b) + plot_layout(widths = c(1.2, 1))) / p6c +
  plot_layout(heights = c(1, 0.9)) +
  plot_annotation(tag_levels = "A")
save_jtm(fig6, "Figure6", h = 150)
saveRDS(list(h6 = h6, lr = lr, nd = nd, od = od),
        "results/figures/final/fig6_data.rds")
# 465. Fig 6A 图例移到下方 ==============================
p6a = p6a +
  guides(fill = guide_colourbar(title = "Spearman ρ with STING score", title.position = "top",
                                direction = "horizontal",
                                barwidth = unit(30, "mm"), barheight = unit(2.5, "mm"))) +
  theme(legend.position = "bottom")

# 466. Fig 6C 左格横轴延长，给文字留空间 ==============================
blank = data.table(test = factor(test_lab["STING_meta"], levels = test_lab), lp = 5.5)
p6c   = p6c + geom_blank(data = blank)

# 467. 重新拼合并导出 Fig 6 ==============================
fig6 = ((p6a | p6b) + plot_layout(widths = c(1, 1.2))) / p6c +
  plot_layout(heights = c(1, 0.9)) +
  plot_annotation(tag_levels = "A")
save_jtm(fig6, "Figure6", h = 150)
saveRDS(list(h6 = h6, lr = lr, nd = nd, od = od),
        "results/figures/final/fig6_data.rds")
# 468. 查看 Fig 1B 与 Fig 7 所需数据 ==============================
file.exists("data_proc/fig1b_data.rds")
mr13 = list.files("results", pattern = "mr_.*\\.csv$", full.names = TRUE)
print(mr13)
trials = fread("results/drug/clinicaltrials_osteosarcoma.csv")
print(colnames(trials))
print(table(trials$relevant, useNA = "ifany"))
print(trials[relevant != "Yes" & relevant != "No", 1:4])
# 469. 用 ClinicalTrials.gov API 核实 4 个待定试验 ==============================
library(jsonlite)
chk_ids = c("NCT04383288", "NCT01194284", "NCT01459484", "NCT02517918")
get_trial = function(id) {
  f = tempfile(fileext = ".json")
  system2("curl.exe", c("-s", "--ssl-no-revoke", "-o", f,
                        paste0("https://clinicaltrials.gov/api/v2/studies/", id)))
  j  = fromJSON(f, simplifyVector = TRUE)
  ps = j$protocolSection
  data.table(nct        = id,
             type       = ps$designModule$studyType,
             phase      = paste(ps$designModule$phases, collapse = "/"),
             status     = ps$statusModule$overallStatus,
             conditions = paste(ps$conditionsModule$conditions, collapse = "; "),
             interventions = paste(ps$armsInterventionsModule$interventions$name, collapse = "; "))
}
chk = rbindlist(lapply(chk_ids, get_trial), fill = TRUE)
print(chk[, .(nct, type, phase, status)])
for (i in seq_len(nrow(chk))) {
  cat("\n====", chk$nct[i], "====\n")
  cat("Conditions:   ", chk$conditions[i], "\n")
  cat("Interventions:", chk$interventions[i], "\n")
}
# 470. 更新待定试验，并复查所有"Yes"试验的研究类型 ==============================
trials[nct %in% c("NCT01459484", "NCT02517918"), relevant := "Yes"]
trials[nct %in% c("NCT04383288", "NCT01194284"), relevant := "No_observational"]
yes_ids = trials[grepl("^Yes", relevant), unique(nct)]
yes_chk = rbindlist(lapply(yes_ids, get_trial), fill = TRUE)
print(yes_chk[, .(nct, type, phase)])
obs_ids = yes_chk[type != "INTERVENTIONAL", nct]
trials[nct %in% obs_ids, relevant := "No_observational"]
fwrite(trials, "results/drug/clinicaltrials_osteosarcoma_verified.csv")

# 471. 统计各类药物的骨肉瘤干预性试验数 ==============================
drug_lab = c(mifamurtide = "Mifamurtide (NOD2; macrophage activation)",
             zoledronic  = "Zoledronic acid (osteoclast)",
             denosumab   = "Denosumab (RANKL; osteoclast)",
             pexidartinib = "CSF1R inhibitor (surufatinib)",
             manganese   = "Manganese",
             sting       = "STING agonists",
             cilengitide = "Integrin αvβ3 inhibitor")
cnt = trials[, .(n = uniqueN(nct[grepl("^Yes", relevant)])), by = query]
cnt[, drug := drug_lab[query]]
print(cnt)
# 472. 重画 Fig 7B：骨肉瘤干预性试验数（按靶细胞着色） ==============================
cnt[, compartment := fcase(query %in% c("mifamurtide", "pexidartinib"), "Myeloid",
                           query %in% c("zoledronic", "denosumab", "cilengitide"), "Osteoclast",
                           query %in% c("sting", "manganese"), "Endothelial")]
cnt[, drug := factor(drug, levels = cnt[order(n, decreasing = FALSE), drug])]
p7b = ggplot(cnt, aes(n, drug, fill = compartment)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = n), hjust = -0.4, size = 8 / .pt) +
  scale_fill_manual(values = c(Myeloid = "#009E73", Osteoclast = "#CC79A7", Endothelial = "#D55E00"),
                    name = "Target compartment") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15)), breaks = 0:5) +
  labs(x = "Interventional osteosarcoma trials (ClinicalTrials.gov)", y = NULL) +
  theme(legend.position = "bottom")

# 473. 重画 Fig 7A：可成药靶点在各细胞中的表达 ==============================
sc    = readRDS("data_proc/sc_tme_scored.rds")
ct_order = c("Myeloid", "Endothelial", "Pericyte", "Osteoclast", "Fibroblast/MSC",
             "T/NK", "Chondroblastic OS", "Osteoblastic OS")
sc$ct = factor(unname(ct_lab[as.character(sc$celltype)]), levels = rev(ct_order))
tg = list("Myeloid-directed"    = c("NOD2", "CSF1R", "CD74", "MIF"),
          "Osteoclast-directed" = c("TNFRSF11A", "FDPS", "ITGAV", "ITGB3"),
          "cGAS–STING"          = c("CGAS", "TMEM173", "STING1", "TBK1", "IRF3"),
          "Mn enzymes"          = c("SOD2", "GLUL", "ARG2"))
tg = lapply(tg, intersect, rownames(sc))
print(tg)
p7a = Seurat::DotPlot(sc, features = tg, group.by = "ct",
                      cols = c("grey85", "#3B75AF"), dot.scale = 3.5) +
  theme_jtm() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic")) +
  guides(colour = guide_colourbar(title = "Scaled mean\nexpression"),
         size   = guide_legend(title = "% expressing")) +
  labs(x = NULL, y = NULL)

# 474. 拼合并导出 Fig 7 ==============================
fig7 = p7a / p7b +
  plot_layout(heights = c(1.2, 0.8)) +
  plot_annotation(tag_levels = "A")
save_jtm(fig7, "Figure7", h = 150)
saveRDS(list(cnt = cnt, tg = tg), "results/figures/final/fig7_data.rds")
rm(sc); gc()
# 475. 修正 Fig 7A 分组标题与 STING 基因名 ==============================
p7a$data$feature.groups = factor(p7a$data$feature.groups, levels = names(tg),
                                 labels = c("Myeloid", "Osteoclast", "cGAS–STING", "Mn enzymes"))
p7a = p7a + scale_x_discrete(labels = function(x) sub("TMEM173", "STING1", x))

# 476. 重新拼合并导出 Fig 7，然后释放内存 ==============================
fig7 = p7a / p7b +
  plot_layout(heights = c(1.2, 0.8)) +
  plot_annotation(tag_levels = "A")
save_jtm(fig7, "Figure7", h = 150)
rm(p7a, fig7); gc()
# 477. 读入 R13 的 MR 数据并查找锰特异 SNP 的定义 ==============================
library(TwoSampleMR)
mr13 = readRDS("data_proc/mr_dat_finngen_r13.rds")
str(mr13, max.level = 1)
print(sapply(mr13, function(d) c(n_snp = nrow(d), n_keep = sum(d$mr_keep))))
main = readLines("scripts/00_main_all_lessons.R", encoding = "UTF-8")
hit  = grep("specific|mn_spec|spec_snp", main, ignore.case = TRUE)
print(hit)
cat(main[hit], sep = "\n")
# 478. 查找 SLC30A10 / SLC39A8 对应的 SNP ==============================
hit2 = grep("SLC30A10|SLC39A8|rs10863510|rs13107325|mn_snp", main)
print(hit2)
for (i in hit2) { cat("---- line", i, "----\n"); cat(main[max(1, i - 2):(i + 4)], sep = "\n") }
print(mr13$os[, c("SNP", "effect_allele.exposure", "beta.exposure", "beta.outcome", "se.outcome")])
# 479. 重跑 R13 的 MR 并保存结果 ==============================
mn_spec = c("rs13107325", "rs10863510")
run_mr  = function(d, label) {
  r = generate_odds_ratios(mr(d, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median")))
  s = generate_odds_ratios(mr(d[d$SNP %in% mn_spec, ], method_list = "mr_ivw"))
  s$method = "IVW, Mn-specific SNPs"
  out = as.data.table(rbind(r, s))[, .(method, nsnp, or, or_lci95, or_uci95, pval)]
  out[, outcome := label]
  out
}
mr13_res = rbind(run_mr(mr13$os,   "Osteosarcoma"),
                 run_mr(mr13$bone, "Bone cancer"))
mr13_res[, method := sub("Inverse variance weighted", "IVW", method)]
mr13_res[, method := sub("MR Egger", "MR-Egger", method)]
print(mr13_res[, .(outcome, method, nsnp, or = round(or, 2), lci = round(or_lci95, 2),
                   uci = round(or_uci95, 2), p = signif(pval, 2))])
fwrite(mr13_res, "results/mr_finngen_r13.csv")
# 480. 重画 Fig 1B：MR 森林图（按 JTM 规格） ==============================
fb = mr13_res[outcome == "Osteosarcoma" | method == "IVW"]
fb[, grp := fifelse(outcome == "Osteosarcoma", "Osteosarcoma\n(134 cases)", "Bone cancer\n(338 cases)")]
fb[, row := paste0(method, " (", nsnp, " SNPs)")]
fb[, row := factor(row, levels = rev(unique(row)))]
fb[, txt := sprintf("%.2f (%.2f–%.2f)   %.2f", or, or_lci95, or_uci95, pval)]
fb[, grp := factor(grp, levels = c("Osteosarcoma\n(134 cases)", "Bone cancer\n(338 cases)"))]
p1b = ggplot(fb, aes(or, row)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.3, colour = "grey40") +
  geom_errorbarh(aes(xmin = or_lci95, xmax = or_uci95), height = 0, linewidth = 0.4) +
  geom_point(size = 1.6, shape = 15) +
  geom_text(aes(x = Inf, label = txt), hjust = -0.08, size = 8 / .pt) +
  facet_grid(grp ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4)) +
  coord_cartesian(xlim = c(0.25, 4), clip = "off") +
  labs(x = "OR per SD increase in genetically predicted blood Mn (log scale)", y = NULL,
       title = "Outcome: FinnGen R13                                                         OR (95% CI)        P") +
  theme(strip.placement   = "outside",
        strip.text.y.left = element_text(angle = 0, face = "bold"),
        plot.title        = element_text(size = 8, face = "bold", hjust = 1),
        plot.margin       = margin(5, 85, 5, 5))

# 481. 读入研究设计图并拼合 Fig 1 ==============================
library(magick)
img1a = image_read_pdf("results/figures/Fig1A_study_design.pdf", density = 300)
p1a   = wrap_elements(full = grid::rasterGrob(as.raster(img1a), interpolate = TRUE))
fig1  = p1a / p1b +
  plot_layout(heights = c(121.6, 55)) +
  plot_annotation(tag_levels = "A")
save_jtm(fig1, "Figure1", h = 190)
# 483. 修正 Fig 1B：列标题与数字对齐 ==============================
fb[, or_txt := sprintf("%.2f (%.2f–%.2f)", or, or_lci95, or_uci95)]
fb[, p_txt  := sprintf("%.2f", pval)]
hdr = data.table(grp = factor(levels(fb$grp)[1], levels = levels(fb$grp)),
                 row = factor(levels(fb$row)[length(levels(fb$row))], levels = levels(fb$row)))
x_or = 5.5; x_p = 17
p1b = ggplot(fb, aes(or, row)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.3, colour = "grey40") +
  geom_errorbar(aes(xmin = or_lci95, xmax = or_uci95), width = 0, linewidth = 0.4, orientation = "y") +
  geom_point(size = 1.6, shape = 15) +
  geom_text(aes(x = x_or, label = or_txt), hjust = 0, size = 8 / .pt) +
  geom_text(aes(x = x_p,  label = p_txt),  hjust = 0, size = 8 / .pt) +
  geom_text(data = hdr, aes(x = x_or, y = row, label = "OR (95% CI)"), hjust = 0,
            nudge_y = 0.8, size = 8 / .pt, fontface = "bold") +
  geom_text(data = hdr, aes(x = x_p,  y = row, label = "P"), hjust = 0,
            nudge_y = 0.8, size = 8 / .pt, fontface = "bold") +
  facet_grid(grp ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4)) +
  coord_cartesian(xlim = c(0.25, 4), clip = "off") +
  labs(x = "OR per SD increase in genetically predicted blood Mn (log scale)", y = NULL) +
  theme(strip.placement   = "outside",
        strip.text.y.left = element_text(angle = 0, face = "bold"),
        plot.margin       = margin(15, 85, 5, 5))

fig1 = p1a / p1b +
  plot_layout(heights = c(121.6, 55)) +
  plot_annotation(tag_levels = "A")
save_jtm(fig1, "Figure1", h = 190)
saveRDS(list(mr13_res = mr13_res, fb = fb), "results/figures/final/fig1_data.rds")
# 484. 修正 Fig 1B：P 值列位置与右边距 ==============================
x_p = 12
p1b = p1b +
  theme(plot.margin = margin(15, 140, 5, 5))
p1b$layers[[5]]$aes_params = NULL
p1b$layers[[5]]$mapping$x  = quote(x_p)
p1b$layers[[7]]$mapping$x  = quote(x_p)
fig1 = p1a / p1b +
  plot_layout(heights = c(121.6, 55)) +
  plot_annotation(tag_levels = "A")
save_jtm(fig1, "Figure1", h = 190)
# 484. 修正 Fig 1B（修正版）：重新生成整张图 ==============================
x_or = 5.5; x_p = 12
p1b = ggplot(fb, aes(or, row)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.3, colour = "grey40") +
  geom_errorbar(aes(xmin = or_lci95, xmax = or_uci95), width = 0, linewidth = 0.4, orientation = "y") +
  geom_point(size = 1.6, shape = 15) +
  geom_text(aes(x = x_or, label = or_txt), hjust = 0, size = 8 / .pt) +
  geom_text(aes(x = x_p,  label = p_txt),  hjust = 0, size = 8 / .pt) +
  geom_text(data = hdr, aes(x = x_or, y = row, label = "OR (95% CI)"), hjust = 0,
            nudge_y = 0.8, size = 8 / .pt, fontface = "bold") +
  geom_text(data = hdr, aes(x = x_p,  y = row, label = "P"), hjust = 0,
            nudge_y = 0.8, size = 8 / .pt, fontface = "bold") +
  facet_grid(grp ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4)) +
  coord_cartesian(xlim = c(0.25, 4), clip = "off") +
  labs(x = "OR per SD increase in genetically predicted blood Mn (log scale)", y = NULL) +
  theme(strip.placement   = "outside",
        strip.text.y.left = element_text(angle = 0, face = "bold"),
        plot.margin       = margin(15, 140, 5, 5))

fig1 = p1a / p1b +
  plot_layout(heights = c(121.6, 55)) +
  plot_annotation(tag_levels = "A")
save_jtm(fig1, "Figure1", h = 190)
# 485. Fig 1B 拆成"森林图 + 数字表"两部分 ==============================
p1b_plot = ggplot(fb, aes(or, row)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.3, colour = "grey40") +
  geom_errorbar(aes(xmin = or_lci95, xmax = or_uci95), width = 0, linewidth = 0.4, orientation = "y") +
  geom_point(size = 1.6, shape = 15) +
  facet_grid(grp ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4)) +
  coord_cartesian(xlim = c(0.25, 4)) +
  labs(x = "OR per SD increase in genetically\npredicted blood Mn (log scale)", y = NULL) +
  theme(strip.placement   = "outside",
        strip.text.y.left = element_text(angle = 0, face = "bold"))

p1b_tab = ggplot(fb, aes(y = row)) +
  geom_text(aes(x = 0,    label = or_txt), hjust = 0, size = 8 / .pt) +
  geom_text(aes(x = 0.78, label = p_txt),  hjust = 0, size = 8 / .pt) +
  geom_text(data = hdr, aes(x = 0,    y = row, label = "OR (95% CI)"), hjust = 0,
            nudge_y = 0.8, size = 8 / .pt, fontface = "bold") +
  geom_text(data = hdr, aes(x = 0.78, y = row, label = "P"), hjust = 0,
            nudge_y = 0.8, size = 8 / .pt, fontface = "bold") +
  facet_grid(grp ~ ., scales = "free_y", space = "free_y") +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme_void(base_size = 8, base_family = "Arial") +
  theme(strip.text = element_blank(), plot.margin = margin(15, 5, 5, 8))

p1b  = (p1b_plot | p1b_tab) + plot_layout(widths = c(2.3, 1))
fig1 = p1a / p1b +
  plot_layout(heights = c(121.6, 55)) +
  plot_annotation(tag_levels = "A")
save_jtm(fig1, "Figure1", h = 190)
# 486. 修正 Fig 1B：列标题不占纵轴，去掉多余标号 ==============================
hdr2 = data.table(grp = factor(levels(fb$grp)[1], levels = levels(fb$grp)))
p1b_tab = ggplot(fb, aes(y = row)) +
  geom_text(aes(x = 0,    label = or_txt), hjust = 0, size = 8 / .pt) +
  geom_text(aes(x = 0.78, label = p_txt),  hjust = 0, size = 8 / .pt) +
  geom_text(data = hdr2, aes(x = 0,    y = Inf, label = "OR (95% CI)"), hjust = 0,
            vjust = -0.5, size = 8 / .pt, fontface = "bold", inherit.aes = FALSE) +
  geom_text(data = hdr2, aes(x = 0.78, y = Inf, label = "P"), hjust = 0,
            vjust = -0.5, size = 8 / .pt, fontface = "bold", inherit.aes = FALSE) +
  facet_grid(grp ~ ., scales = "free_y", space = "free_y") +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme_void(base_size = 8, base_family = "Arial") +
  theme(strip.text = element_blank(), plot.margin = margin(15, 5, 5, 8))

p1b  = wrap_elements(full = (p1b_plot | p1b_tab) + plot_layout(widths = c(2.3, 1)))
fig1 = p1a / p1b +
  plot_layout(heights = c(121.6, 55)) +
  plot_annotation(tag_levels = "A")
save_jtm(fig1, "Figure1", h = 190)
saveRDS(list(mr13_res = mr13_res, fb = fb), "results/figures/final/fig1_data.rds")
# 487. 加入 TMEM173 重算 T4 评分 ==============================
library(Seurat)
sc     = readRDS("data_proc/sc_tme_scored.rds")
g4b    = intersect(c("CGAS", "TMEM173", "TBK1", "IRF3", "IFNB1", "CXCL10", "ISG15"), rownames(sc))
print(g4b)
sc     = AddModuleScore(sc, features = list(g4b), name = "T4b")
md     = as.data.table(sc@meta.data)[, .(patient, celltype, T4 = Tier4, T4b = T4b1)]
print(round(cor(md$T4, md$T4b), 3))                       # 新旧评分的相关性

# 488. 伪批量配对比较：原版 vs 加入 STING1 ==============================
pb4  = md[, .(n = .N, T4 = mean(T4), T4b = mean(T4b)), by = .(patient, celltype)][n >= 20]
cmp  = function(cell, var) {
  w = dcast(pb4, patient ~ celltype, value.var = var)
  d = na.omit(w[[cell]] - w[["OS_Osteoblastic"]])
  data.table(cell = cell, score = var, k = sum(d > 0), n = length(d),
             p_sign = binom.test(sum(d > 0), length(d))$p.value)
}
sens = rbindlist(lapply(c("Endothelial", "Pericyte", "T_NK", "Myeloid"), function(cl)
  rbind(cmp(cl, "T4"), cmp(cl, "T4b"))))
print(sens)
fwrite(sens, "results/sensitivity_sting1_discovery.csv")
rm(sc); gc()
# 489. 从主脚本中查找 Methods 所需参数 ==============================
main = readLines("scripts/00_main_all_lessons.R", encoding = "UTF-8")
show = function(pat, n = 2) {
  for (i in grep(pat, main)) { cat("---- line", i, "----\n"); cat(main[i:(i + n)], sep = "\n") }
}
show("ieugwas|gwasinfo|extract_instruments|clump", 3)  # 锰 GWAS 编号、P 阈值、clump 参数
show("mad|MAD|nFeature|percent.mt", 2)                 # QC 阈值
show("FindVariableFeatures|RunPCA|dims|resolution", 1) # 可变基因数、PC 数、分辨率
show("computeCommunProb|filterCommunication|CellChatDB", 1)
show("tumor_rich|tumour_rich|quantile", 2)             # 空间肿瘤富集区的定义
show("coxph", 1)                                       # Cox 模型的协变量
show("ARG1", 0)                                        # T2 是否含 ARG1
packageVersion("Seurat"); packageVersion("TwoSampleMR"); packageVersion("CellChat")
# 490. 确认第 166 步的基因集修改，以及工具变量的 P 阈值和来源 ==============================
cat(main[752:768], sep = "\n")
show("5e-8|5e-08|1e-5|1e-05|P <|P<", 1)
show("mn_gwas|download|zenodo|gwas_catalog|GCST", 1)
# 491. 确认 Tier 评分的计算位置，以及锰 GWAS 的来源 ==============================
show('name = "Tier"', 1)
cat(main[40:68], sep = "\n")
# 492. 重画 Fig 2C（补上 STING1）并重新拼合 Fig 2 ==============================
library(Seurat)
source("scripts/15_figure_theme.R", encoding = "UTF-8")
sc    = readRDS("data_proc/sc_tme_scored.rds")
sc$ct = factor(unname(ct_lab[as.character(sc$celltype)]), levels = names(ct_col))
f2    = readRDS("results/figures/final/fig2_data.rds")
mn_genes = list(
  "T1 transport"     = c("SLC39A8", "SLC39A14", "SLC30A10", "SLC11A2", "ATP2C1", "ATP13A1", "TMEM165", "TFRC"),
  "T2 Mn enzymes"    = c("SOD2", "ARG2", "GLUL", "PC", "PPM1A", "PPM1B"),
  "T3 glycosylation" = c("B4GALT1", "GALNT1", "GALNT2", "XYLT1", "B3GAT3", "MGAT1"),
  "T4 cGAS–STING"    = c("CGAS", "TMEM173", "TBK1", "IRF3", "IFNB1", "CXCL10", "ISG15"))
p2a = plot_umap(sc)
p2b = ggplot(f2$hm, aes(tier, ct, fill = z)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", z)), size = 8 / .pt) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, name = "z-score\n(within tier)") +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.line = element_blank(), axis.ticks = element_blank())
p2c = DotPlot(sc, features = mn_genes, group.by = "ct", cols = c("grey85", "#3B75AF"), dot.scale = 3.5) +
  theme_jtm() +
  scale_x_discrete(labels = function(x) sub("TMEM173", "STING1", x)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
        legend.position = "bottom", legend.box = "horizontal", legend.margin = margin(0, 0, 0, 0)) +
  guides(colour = guide_colourbar(title = "Scaled mean expression", title.position = "top",
                                  barwidth = unit(25, "mm"), barheight = unit(2.5, "mm")),
         size   = guide_legend(title = "% expressing", title.position = "top")) +
  labs(x = NULL, y = NULL)
les_col = c("Primary" = "#0072B2", "Primary or recurrent" = "#56B4E9",
            "Recurrent" = "#E69F00", "Lung metastasis" = "#CC79A7")
p2d = ggplot(f2$pb, aes(score, ct)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.3, width = 0.6, fill = NA) +
  geom_point(aes(colour = lesion), size = 0.8, position = position_jitter(height = 0.15, seed = 1)) +
  facet_wrap(~ tier, nrow = 1, scales = "free_x") +
  scale_colour_manual(values = les_col, name = "Lesion") +
  scale_x_continuous(breaks = scales::breaks_pretty(n = 3)) +
  labs(x = "Mean module score per patient", y = NULL) +
  theme(legend.position = "bottom", panel.grid.major.x = element_line(linewidth = 0.2, colour = "grey90"))
fig2 = (free(p2a) | p2b) / p2c / p2d +
  plot_layout(heights = c(1.1, 1, 1)) +
  plot_annotation(tag_levels = "A")
save_jtm(fig2, "Figure2", h = 220)
rm(sc, p2a, p2c, fig2); gc()
# 493. 统一 Fig 2C 行顺序，加大 Fig 2D 分格间距 ==============================
ct_order = c("Myeloid", "Endothelial", "Pericyte", "Osteoclast", "Fibroblast/MSC",
             "T/NK", "Chondroblastic OS", "Osteoblastic OS")
sc    = readRDS("data_proc/sc_tme_scored.rds")
sc$ct = factor(unname(ct_lab[as.character(sc$celltype)]), levels = rev(ct_order))
p2a   = plot_umap(sc)
p2c   = DotPlot(sc, features = mn_genes, group.by = "ct", cols = c("grey85", "#3B75AF"), dot.scale = 3.5) +
  theme_jtm() +
  scale_x_discrete(labels = function(x) sub("TMEM173", "STING1", x)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
        legend.position = "bottom", legend.box = "horizontal", legend.margin = margin(0, 0, 0, 0)) +
  guides(colour = guide_colourbar(title = "Scaled mean expression", title.position = "top",
                                  barwidth = unit(25, "mm"), barheight = unit(2.5, "mm")),
         size   = guide_legend(title = "% expressing", title.position = "top")) +
  labs(x = NULL, y = NULL)
p2d  = p2d + theme(panel.spacing.x = unit(4, "mm"))
fig2 = (free(p2a) | p2b) / p2c / p2d +
  plot_layout(heights = c(1.1, 1, 1)) +
  plot_annotation(tag_levels = "A")
save_jtm(fig2, "Figure2", h = 220)
rm(sc, p2a, p2c, fig2); gc()
# 494. 提取 Methods 中剩余的参数 ==============================
iv = fread("data_proc/mn_iv_main_clumped.csv")
print(iv[, .(rsid, F)]); print(round(c(mean_F = mean(iv$F), min_F = min(iv$F)), 1))
show("collapse|avereps|duplicated|lumiHumanAll", 2)          # GSE21257 重复探针如何合并
show("finngen_R13_C3_OSTEOSARCOMA|finngen_R5|ncontrol|n_control", 1)
pk = c("scDblFinder", "harmony", "MRPRESSO", "MVMR", "survival", "ggplot2",
       "patchwork", "lumiHumanAll.db", "ieugwasr")
print(sapply(pk, function(p) tryCatch(as.character(packageVersion(p)), error = function(e) NA)))
# 495. 提取 Results 所需的剩余数字 ==============================
print(fread("results/mvmr_bone_iron.csv"))                       # MVMR 结果
print(fread("results/random_geneset_gse21257.csv"))              # 观察 P 与经验 P
print(fread("results/scRNA/discovery_validation_combined.csv")[,
                                                                        .(cell, tier, k_all, n_all, p_sign, p_sign_BH)])                     # 软骨母细胞型 P 是否 BH 校正
show("p_meta =|p_meta=|wilcox.test\\(.*meta", 1)                           # STING-转移用的检验方法
g21 = fread("data_proc/gse21257_clin_scores.csv")
library(survival)
print(summary(coxph(Surv(os_mon, OS) ~ scale(T3_Glycosylation), data = g21))$conf.int)
print(summary(coxph(Surv(os_mon, OS) ~ scale(T3_Glycosylation) + huvos_good + meta_dx + age_y,
                    data = g21))$conf.int[1, ])
# 496. 打印随机基因集分析的原始代码 ==============================
cat(main[2268:2300], sep = "\n")
# 497. 加载包 =================================================
library(ieugwasr)
library(data.table)

# 498. 检查令牌是否已设置（只显示 TRUE/FALSE，不显示内容） ====
nchar(Sys.getenv("OPENGWAS_JWT")) > 0

# 499. 读取 FinnGen R13 清单文件 ===============================
url13 = "https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_manifest.tsv"
man13 = fread(url13)
names(man13)

# 500. 查骨肉瘤和骨癌的病例数与对照数 ==========================
man13[phenocode %in% c("C3_OSTEOSARCOMA_WIDE", "C3_BONE_CARTILAGE"),
      .(phenocode, num_cases, num_controls)]

# 501. 查 R5 骨癌和铁 GWAS 的信息 ==============================
info = gwasinfo(c("finn-b-C3_BONE_CARTILAGE_EXALLC", "ieu-a-1049"))
keep = intersect(c("id", "trait", "author", "year", "pmid",
                   "sample_size", "ncase", "ncontrol"), names(info))
print(as.data.frame(info)[, keep])

# 502. 计算代理 SNP 的 r² ======================================
ld = ld_matrix(c("rs10863510", "rs6541114"), pop = "EUR")
round(ld^2, 3)

# 503. 定位 MVMR 和生存分析的代码行 ============================
main = readLines("scripts/00_main_all_lessons.R", encoding = "UTF-8")
idx  = grep("mvmr|mv_harmon|strength_|coxph|Surv\\(", main, ignore.case = TRUE)
cat(paste(idx, main[idx]), sep = "\n")

# 504. 查药物检索结果文件的保存日期 ============================
fi = file.info(list.files("results/drug", full.names = TRUE))
fi[, "mtime", drop = FALSE]
# 506. 检查 MVMR 对象是否还在内存中 ============================
exists("mvdat")
exists("mv_in2")

# 507. 如果都是 TRUE，直接计算条件 F 统计量 ====================
library(MVMR)
strength_mvmr(mv_in2, gencov = 0)

# 508. 查看 TARGET 生存分析代码的上下文 ========================
cat(paste(1170:1200, main[1170:1200]), sep = "\n")
# 509. 打印 MVMR 那一段代码（第 351–474 行） ==================
cat(paste(351:474, main[351:474]), sep = "\n")

# 510. 打印 TARGET 生存分析的上下文（第 1170–1200 行） ========
cat(paste(1170:1200, main[1170:1200]), sep = "\n")
# 511. 查 iv、dat、bs、cl2、bulk_sets 在哪几行生成或保存 =======
idx = grep("^(iv|dat|bs|cl2|bulk_sets|bk)\\s*=|saveRDS\\(|fwrite\\(", main)
cat(paste(idx, main[idx]), sep = "\n")

# 512. 列出已经保存的 rds 和 csv 文件 ==========================
list.files("data_proc", pattern = "\\.(rds|csv)$")
list.files("results",   pattern = "\\.(rds|csv)$", recursive = TRUE)

# 513. 查看 Fig. 6 数据文件里保存了什么 ========================
f6 = readRDS("results/figures/final/fig6_data.rds")
names(f6)
# 514. 设置工作目录并加载包 ====================================
# [local-only, disabled for public release] setwd(".")
library(TwoSampleMR)
library(ieugwasr)
library(data.table)
library(MVMR)

# 515. 读取保存好的主 MR 数据和锰工具变量 ======================
dat = readRDS("data_proc/mr_dat_bone.rds")
iv  = fread("data_proc/mn_iv_main_clumped.csv")

# 516. 提取血清铁的工具变量 ====================================
fe_iv = extract_instruments("ieu-a-1049")

# 517. 合并锰和铁的 SNP，取较小的 P 值，重新 clumping ===========
mv_in    = rbind(data.frame(rsid = dat$SNP,   pval = dat$pval.exposure),
                 data.frame(rsid = fe_iv$SNP, pval = fe_iv$pval.exposure))
mv_in    = aggregate(pval ~ rsid, data = mv_in, FUN = min)
mv_in$id = "mvmr"
mv_clump = ld_clump(mv_in, clump_kb = 10000, clump_r2 = 0.001, pop = "EUR")
snps_mv  = mv_clump$rsid
length(snps_mv)

# 518. 从完整锰 GWAS 中取出这些 SNP =============================
pos_tab = rbind(data.frame(SNP = iv$rsid, chrpos = paste0(iv$CHR, ":", iv$POS)),
                data.frame(SNP = "rs1525892", chrpos = "3:133484712"))
mn = readRDS("data_proc/mn_gwas_clean.rds")
setDT(mn)
mn[, chrpos := paste0(CHR, ":", POS)]
mn_mv = merge(mn[chrpos %in% pos_tab$chrpos], pos_tab, by = "chrpos")
rm(mn); gc()

# 519. 整理成锰的暴露格式 ======================================
exp_mn = format_data(as.data.frame(mn_mv), type = "exposure",
                     snp_col = "SNP", beta_col = "Effect", se_col = "StdErr",
                     effect_allele_col = "Allele1", other_allele_col = "Allele2",
                     eaf_col = "Freq1", pval_col = "P")
exp_mn$exposure = "Manganese"; exp_mn$id.exposure = "Mn"

# 520. 提取铁的效应和骨癌结局（都允许代理 SNP） ================
fe_out = extract_outcome_data(snps = snps_mv, outcomes = "ieu-a-1049", proxies = TRUE)
exp_fe = convert_outcome_to_exposure(fe_out)
exp_fe$exposure = "Iron"; exp_fe$id.exposure = "Fe"
out_mv = extract_outcome_data(snps = snps_mv,
                              outcomes = "finn-b-C3_BONE_CARTILAGE_EXALLC", proxies = TRUE)

# 521. 统一方向并重跑 MVMR，确认和原来的结果一致 ===============
mv_exp = as.data.frame(rbindlist(list(exp_mn, exp_fe), fill = TRUE))
mvdat  = mv_harmonise_data(mv_exp, out_mv)
mv_res = mv_multiple(mvdat)
mv_res$result[, c("exposure", "nsnp", "b", "se", "pval")]
fread("results/mvmr_bone_iron.csv")

# 522. 计算条件 F 统计量和异质性 Q =============================
mv_in2 = format_mvmr(BXGs   = mvdat$exposure_beta,
                     BYG    = mvdat$outcome_beta,
                     seBXGs = mvdat$exposure_se,
                     seBYG  = mvdat$outcome_se,
                     RSID   = rownames(mvdat$exposure_beta))
strength_mvmr(mv_in2, gencov = 0)
pleiotropy_mvmr(mv_in2, gencov = 0)
# 523. 读取 TARGET 评分和临床数据，确认列名 =====================
bk = fread("data_proc/target_clin_scores.csv")
names(bk)

# 524. 找出四个锰评分的列名 ====================================
tiers = grep("^T[0-9]", names(bk), value = TRUE)
tiers

# 525. 单因素 Cox 回归（和原来第 1187–1196 行相同） ============
library(survival)
cox_one = function(tier, time, event) {
  d = bk[!is.na(get(time)) & !is.na(get(event))]
  f = as.formula(paste0("Surv(", time, ", ", event, ") ~ scale(", tier, ")"))
  s = summary(coxph(f, data = d))
  data.table(tier, endpoint = event, n = nrow(d), events = sum(d[[event]]),
             HR = round(s$conf.int[1], 2), lo = round(s$conf.int[3], 2),
             hi = round(s$conf.int[4], 2), p = signif(s$coefficients[5], 3))
}
tg_cox = rbindlist(c(lapply(tiers, cox_one, time = "os_days",  event = "OS"),
                     lapply(tiers, cox_one, time = "efs_days", event = "EFS")))
tg_cox

# 526. 保存结果 ================================================
fwrite(tg_cox, "results/target_cox_tiers.csv")
# 527. 确认 exposure1 和 exposure2 分别是锰还是铁 ==============
colnames(mvdat$exposure_beta)
# 528. 加载作图主题和包 ========================================
# [local-only, disabled for public release] setwd(".")
source("scripts/15_figure_theme.R", encoding = "UTF-8")
library(survival)
library(ggplot2)
library(data.table)

# 529. 读取 TARGET 结果和 GSE21257 数据 ========================
tg  = fread("results/target_cox_tiers.csv")
g21 = fread("data_proc/gse21257_clin_scores.csv")
tiers = c("T1_Transport", "T2_Enzyme", "T3_Glycosylation", "T4_STING")

# 530. GSE21257 单因素 Cox 回归（总生存） ======================
cox_g = function(tier) {
  f = as.formula(paste0("Surv(os_mon, OS) ~ scale(", tier, ")"))
  s = summary(coxph(f, data = g21))
  data.table(tier, endpoint = "OS", n = s$n, events = s$nevent,
             HR = s$conf.int[1], lo = s$conf.int[3], hi = s$conf.int[4],
             p  = s$coefficients[5])
}
g_cox = rbindlist(lapply(tiers, cox_g))
g_cox

# 531. 合并两个队列并加标签 ====================================
tg[, cohort := "TARGET-OS"]
g_cox[, cohort := "GSE21257"]
fs5 = rbind(tg, g_cox, fill = TRUE)
nm  = c(T1_Transport = "Mn transport", T2_Enzyme = "Mn-dependent enzymes",
        T3_Glycosylation = "Mn glycosylation", T4_STING = "cGAS\u2013STING")
fs5[, tier_name := factor(nm[tier], levels = rev(nm))]
fs5[, ep_name   := ifelse(endpoint == "OS", "overall survival", "event-free survival")]
fs5[, panel     := sprintf("%s: %s\n(n = %d, %d events)", cohort, ep_name, n, events)]
fs5[, panel     := factor(panel, levels = unique(panel[order(cohort != "TARGET-OS", endpoint != "OS")]))]
fs5[, sig       := p < 0.05]
fs5[, lbl       := sprintf("%.2f (%.2f\u2013%.2f), P = %s", HR, lo, hi,
                           ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))]
fs5[, .(panel, tier_name, lbl)]

# 532. 画森林图 ================================================
p_s5 = ggplot(fs5, aes(x = HR, y = tier_name)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.4, colour = "grey50") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2, linewidth = 0.5) +
  geom_point(aes(fill = sig), shape = 21, size = 2, stroke = 0.5) +
  geom_text(aes(x = 3.2, label = lbl), hjust = 0, size = 8 / .pt) +
  scale_fill_manual(values = c(`FALSE` = "white", `TRUE` = "#C0392B"), guide = "none") +
  scale_x_log10(limits = c(0.4, 12), breaks = c(0.5, 1, 2)) +
  facet_wrap(~ panel, nrow = 1) +
  labs(x = "Hazard ratio per SD (log scale)", y = NULL) +
  theme_jtm() +
  theme(strip.background = element_blank(),
        strip.text       = element_text(size = 8, face = "bold"))
p_s5

# 533. 保存图和作图数据 ========================================
save_jtm(p_s5, "FigS5_survival_forest", h = 65)
saveRDS(fs5, "results/figures/final/figS5_data.rds", compress = "xz")
# 534. 修改森林图：数字放在区间上方，收窄 x 轴 =================
p_s5 = ggplot(fs5, aes(x = HR, y = tier_name)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.4, colour = "grey50") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.15, linewidth = 0.5) +
  geom_point(aes(fill = sig), shape = 21, size = 2, stroke = 0.5) +
  geom_text(aes(x = 1.06, label = lbl), nudge_y = 0.33, size = 8 / .pt) +
  scale_fill_manual(values = c(`FALSE` = "white", `TRUE` = "#C0392B"), guide = "none") +
  scale_x_log10(limits = c(0.4, 2.8), breaks = c(0.5, 1, 2)) +
  scale_y_discrete(expand = expansion(add = c(0.5, 0.8))) +
  facet_wrap(~ panel, nrow = 1) +
  labs(x = "Hazard ratio per SD (log scale)", y = NULL) +
  theme_jtm() +
  theme(strip.background = element_blank(),
        strip.text       = element_text(size = 8, face = "bold"))

# 535. 保存并打开 PDF 检查实际效果 =============================
save_jtm(p_s5, "FigS5_survival_forest", h = 75)
saveRDS(fs5, "results/figures/final/figS5_data.rds", compress = "xz")
# [local-only, disabled for public release] shell.exec("results/figures/final/FigS5_survival_forest.pdf")
# 539. 给数字加白色底，避免虚线穿过文字 ========================
p_s5 = p_s5 + NULL
p_s5$layers[[4]] = NULL
p_s5 = p_s5 +
  geom_label(aes(x = 1.06, label = lbl), nudge_y = 0.33, size = 8 / .pt,
             fill = "white", label.size = 0, label.r = unit(0, "pt"),
             label.padding = unit(0.6, "pt"))

# 540. 保存并打开 PDF ==========================================
save_jtm(p_s5, "FigS5_survival_forest", h = 75)
saveRDS(fs5, "results/figures/final/figS5_data.rds", compress = "xz")
# [local-only, disabled for public release] shell.exec("results/figures/final/FigS5_survival_forest.pdf")
# 541. 查看空间相关结果里的样本名 ==============================
sp = fread("results/spatial_within_sample_cor.csv")
names(sp)
unique(sp[[1]])

# 542. 查看肿瘤区域结果里的样本名和分组 ========================
tp = fread("results/spatial_tumor_spots_by_patient.csv")
names(tp)
head(tp, 30)
# 543. 从样本名中提取真实患者 ID（同一患者编号数字相同） =========
sp[, sample := sub("^GSM[0-9]+_", "", prefix)]
sp[, pid    := paste0("P", gsub("[^0-9]", "", sample))]
sp[, .(n_sections = uniqueN(sample)), by = pid][order(pid)]
uniqueN(sp$pid)

# 544. 每名患者取各切片相关系数的平均，再按配对汇总 ===============
pt_cor = sp[, .(rho = mean(rho)), by = .(pid, cell, mn)]
sum_cor = pt_cor[, .(n_pt     = .N,
                     median   = round(median(rho), 3),
                     n_pos    = sum(rho > 0)), by = .(cell, mn)]
sum_cor[order(cell, mn)]

# 545. 肿瘤区域：加入真实患者 ID ================================
tp[, pid := paste0("P", gsub("[^0-9]", "", patient))]
tp[, .(pid, patient, grp2)][order(pid)]

# 546. 原发 vs 复发：线性混合模型（患者为随机效应） ===============
library(nlme)
tiers_sp = c("T1_Transport", "T2_Enzyme", "T3_Glyco", "T4_STING")
mm = rbindlist(lapply(tiers_sp, function(t) {
  f   = as.formula(paste(t, "~ grp2"))
  fit = lme(f, random = ~ 1 | pid, data = tp, method = "REML")
  tt  = summary(fit)$tTable
  data.table(tier = t, diff = round(tt[2, "Value"], 3), p = signif(tt[2, "p-value"], 3))
}))
mm

# 547. 敏感性分析：5 对配对样本的配对检验 =========================
pairs = tp[pid %in% tp[, .N, by = pid][N == 2, pid]]
pw = rbindlist(lapply(tiers_sp, function(t) {
  w = dcast(pairs, pid ~ grp2, value.var = t)
  data.table(tier = t, n_pairs = nrow(w),
             p = signif(wilcox.test(w$Relapse_Met, w$Primary, paired = TRUE)$p.value, 3))
}))
pw

# 548. 保存重算结果 ============================================
fwrite(sum_cor, "results/spatial_cor_by_true_patient.csv")
fwrite(mm,      "results/spatial_primary_vs_relapse_lme.csv")
# 549. 读取全部"细胞类型 × 锰程序"相关结果，按真实患者重算 ======
sa = fread("results/spatial_all_pairs_cor.csv")
names(sa)
sa[, sample := sub("^GSM[0-9]+_", "", prefix)]
sa[, pid    := paste0("P", gsub("[^0-9]", "", sample))]
uniqueN(sa$pid)
sa_pt  = sa[, .(rho = mean(rho)), by = .(pid, cell, mn)]
sa_sum = sa_pt[, .(median = round(median(rho), 3), n_pos = sum(rho > 0)), by = .(cell, mn)]
sa_sum[grepl("steoblast|umor", cell)]

# 550. 查看 Fig. 5 作图数据的结构 ==============================
f5 = readRDS("results/figures/final/fig5_data.rds")
names(f5)
str(f5, max.level = 1)
# 551. 查看各个数据表的前几行 ==================================
head(f5$hm5, 3); head(f5$pc, 3); f5$pos; head(f5$td, 3); f5$pw

# 552. 找到 Fig. 5 作图代码所在的行 ============================
idx = grep("fig5|p5a|p5b|p5c|p5d|hm5|f5\\$", main)
cat(paste(idx, main[idx]), sep = "\n")
# 553. 在所有脚本里搜索 Fig. 5 的作图代码 =======================
fs = list.files("scripts", pattern = "\\.R$", full.names = TRUE)
for (f in fs) {
  x = readLines(f, encoding = "UTF-8", warn = FALSE)
  i = grep("fig5_data|hm5|pos\\$|pw\\$|td\\$", x)
  if (length(i)) cat("\n=====", basename(f), "=====\n", paste(i, x[i]), sep = "\n")
}

# 554. 查看各个数据表的列名和前几行 ============================
head(f5$hm5, 3); head(f5$pc, 3); f5$pos; head(f5$td, 3); f5$pw
# 555. 细胞类型和锰程序的显示名称（与原图一致） ================
cell_map = c(Myeloid = "Myeloid", Endothelial = "Endothelial", Osteoclast = "Osteoclast",
             Fibroblast = "Fibroblast/MSC", T_NK = "T/NK", OS_Osteoblastic = "Osteoblastic OS")
unique(sa$cell); unique(sa$mn)
levels(f5$hm5$mn); unique(f5$hm5$cell)

# 556. 重建 B 图数据：19 名患者的中位相关系数 ===================
mn_map = setNames(levels(f5$hm5$mn), c("T1_Transport", "T2_Enzyme", "T3_Glyco", "T4_STING", "aVb3"))
hm5_new = sa_pt[cell %in% names(cell_map),
                .(rho = median(rho)), by = .(cell, mn)]
hm5_new[, cell := cell_map[cell]]
hm5_new[, mn   := factor(mn_map[mn], levels = levels(f5$hm5$mn))]
hm5_new[order(cell, mn)]

# 557. 重建 C 图数据：每名患者一个点（19 名） ===================
lab_key = unique(f5$pc[, .(cell, mn, lab)])
pc_new  = merge(pt_cor, lab_key, by = c("cell", "mn"))
setnames(pc_new, "pid", "patient")
pc_new[, .N, by = lab]

# 558. 重建 C 图右侧的阳性人数 ==================================
pos_new = pc_new[, .(k = sum(rho > 0), n = .N, med = median(rho)), by = lab]
pos_new

# 559. 重建 D 图的 P 值（混合模型，另加 BH 校正） ===============
pw_new = copy(f5$pw)
pw_new[, p     := mm$p]
pw_new[, p_adj := p.adjust(mm$p, method = "BH")]
pw_new

# 560. 组合成新的 Fig. 5 数据并保存 ============================
f5_new = list(sc5 = f5$sc5, rep_p = f5$rep_p, hm5 = hm5_new,
              pc = pc_new, pos = pos_new, td = f5$td, pw = pw_new)
saveRDS(f5_new, "results/figures/final/fig5_data_19pt.rds", compress = "xz")
# 561. 加载作图主题和包 ========================================
source("scripts/15_figure_theme.R", encoding = "UTF-8")
library(ggplot2)
library(patchwork)

# 562. B 图：19 名患者的中位相关系数热图 =======================
cell_lv = c("Myeloid", "Endothelial", "Osteoclast", "Fibroblast/MSC", "T/NK", "Osteoblastic OS")
hm = copy(hm5_new)
hm[, cell := factor(cell, levels = rev(cell_lv))]
hm[, txt  := sub("^-0\\.00$", "0.00", sprintf("%.2f", rho))]
p5b = ggplot(hm, aes(mn, cell, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = txt), size = 8 / .pt) +
  scale_fill_gradient2(low = "#4575B4", mid = "white", high = "#B2182B", midpoint = 0,
                       limits = c(-0.15, 0.15), oob = scales::squish,
                       name = "Median Spearman \u03c1") +
  labs(x = NULL, y = NULL) +
  theme_jtm() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line = element_blank(), axis.ticks = element_blank(),
        legend.position = "bottom") +
  guides(fill = guide_colourbar(title.position = "top",
                                barwidth = unit(30, "mm"), barheight = unit(2.5, "mm")))

# 563. C 图：每名患者的相关系数，右侧标出阳性人数 ===============
lv = pos_new[order(-med), as.character(lab)]
pc = copy(pc_new); pos = copy(pos_new)
pc[,  lab := factor(lab, levels = rev(lv))]
pos[, lab := factor(lab, levels = rev(lv))]
set.seed(1)
p5c = ggplot(pc, aes(rho, lab)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  geom_jitter(height = 0.15, width = 0, size = 1.2, colour = "#3C6EA8", alpha = 0.75) +
  geom_errorbar(data = pos, aes(y = lab, xmin = med, xmax = med), inherit.aes = FALSE,
                width = 0.45, colour = "#C0392B", linewidth = 0.6) +
  geom_text(data = pos, aes(x = Inf, y = lab, label = paste0(k, "/", n)),
            inherit.aes = FALSE, hjust = -0.3, size = 8 / .pt) +
  annotate("text", x = Inf, y = Inf, label = "\u03c1 > 0", hjust = -0.3, vjust = -0.6,
           fontface = "bold", size = 8 / .pt) +
  scale_y_discrete(labels = function(x) gsub(" \u2013 ", "\n", x)) +
  coord_cartesian(clip = "off") +
  labs(x = "Spearman \u03c1 per patient", y = NULL) +
  theme_jtm() +
  theme(plot.margin = margin(12, 30, 5, 5, "pt"))

# 564. D 图：原发 vs 复发，配对患者用灰线连接 ===================
td = copy(f5$td)
td[, pid := paste0("P", gsub("[^0-9]", "", patient))]
paired = td[, uniqueN(grp2), by = pid][V1 == 2, pid]
pw = copy(pw_new)
pw[, lab := sprintf("P = %.3f", p)]
p5d = ggplot(td, aes(grp2, score, colour = grp2)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, linewidth = 0.5, fill = NA) +
  geom_line(data = td[pid %in% paired], aes(group = pid), colour = "grey70", linewidth = 0.3) +
  geom_point(size = 1.2) +
  geom_text(data = pw, aes(x = 1.5, y = y * 1.08, label = lab), inherit.aes = FALSE, size = 8 / .pt) +
  facet_wrap(~ tier, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = c("#0072B2", "#CC79A7"), guide = "none") +
  scale_x_discrete(labels = c("Primary", "Relapse/\nmetastasis")) +
  labs(x = NULL, y = "Mean score in\ntumor-rich spots") +
  theme_jtm() +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold", size = 8))

# 565. 预览 B、C、D，并查看 A 图数据的结构 ======================
(p5b | p5c) / p5d
names(f5$sc5); head(f5$sc5, 3); f5$rep_p
# 566. 修正 B 图：图例放右侧，只保留三个刻度 ===================
p5b = p5b +
  scale_fill_gradient2(low = "#4575B4", mid = "white", high = "#B2182B", midpoint = 0,
                       limits = c(-0.15, 0.15), breaks = c(-0.1, 0, 0.1), oob = scales::squish,
                       name = "Median\nSpearman \u03c1") +
  guides(fill = guide_colourbar(barwidth = unit(2.5, "mm"), barheight = unit(22, "mm"))) +
  theme(legend.position = "right")

# 567. 修正 C 图：按正文顺序排列，标签拆成两行 ==================
ord = c("Myeloid", "Osteoclast.*\u03b1v\u03b23", "Endothelial", "Fibroblast", "Osteoclast.*transport", "T/NK")
lv  = sapply(ord, function(o) grep(o, levels(pc_new$lab), value = TRUE))
pc[,  lab := factor(as.character(lab), levels = rev(lv))]
pos[, lab := factor(as.character(lab), levels = rev(lv))]
p5c = p5c %+% pc +
  scale_y_discrete(labels = function(x) gsub("\\s[\u2013\u2014-]\\s", "\n", x))
p5c$layers[[3]]$data = pos
p5c$layers[[4]]$data = pos

# 568. 再次预览，并查看 A 图数据的结构 ==========================
(p5b | p5c) / p5d
names(f5$sc5); head(f5$sc5, 3); f5$rep_p
# 569. A 图：代表性切片的五张空间分布图 =========================
library(ggrastr)
a_vars = c(Tumor = "Tumor (osteoblastic)", Myeloid = "Myeloid", T2 = "Mn enzymes",
           Osteoclast = "Osteoclast", aVb3 = "\u03b1v\u03b23 integrin")
one_map = function(v) {
  d   = f5$sc5[, .(x, y, val = get(v))]
  lim = quantile(d$val, c(0.01, 0.99))
  ggplot(d, aes(x, y, colour = val)) +
    rasterise(geom_point(size = 0.35, stroke = 0), dpi = 300) +
    scale_colour_distiller(palette = "Spectral", limits = lim, oob = scales::squish,
                           breaks = scales::breaks_pretty(3), name = NULL) +
    scale_y_reverse() +
    coord_fixed() +
    labs(title = a_vars[[v]]) +
    theme_void(base_size = 8, base_family = "Arial") +
    theme(plot.title  = element_text(face = "bold", size = 8, hjust = 0.5),
          legend.position = "bottom",
          legend.key.width  = unit(10, "mm"), legend.key.height = unit(1.5, "mm"),
          legend.text = element_text(size = 8))
}
p5a = wrap_plots(lapply(names(a_vars), one_map), nrow = 1)

# 570. 备份旧的 Fig. 5（24 名"患者"版本） =======================
fd = "results/figures/final"
for (ext in c("pdf", "tiff")) {
  file.copy(file.path(fd, paste0("Figure5.", ext)),
            file.path(fd, paste0("Figure5_old24pt.", ext)), overwrite = FALSE)
}

# 571. 拼合四个图并保存 ========================================
fig5 = wrap_elements(full = p5a) / (p5b | p5c) / p5d +
  plot_layout(heights = c(1, 1.25, 1.1)) +
  plot_annotation(tag_levels = "A")
save_jtm(fig5, "Figure5", h = 230)
# [local-only, disabled for public release] shell.exec(file.path(fd, "Figure5.pdf"))
# 572. 修正 A 图：对调 x/y、加大点、缩短色条 =====================
one_map = function(v) {
  d   = f5$sc5[, .(x, y, val = get(v))]
  lim = quantile(d$val, c(0.01, 0.99))
  ggplot(d, aes(y, x, colour = val)) +
    rasterise(geom_point(size = 0.6, stroke = 0), dpi = 300) +
    scale_colour_distiller(palette = "Spectral", limits = lim, oob = scales::squish,
                           breaks = scales::breaks_pretty(2), name = NULL) +
    scale_y_reverse() +
    coord_fixed() +
    labs(title = a_vars[[v]]) +
    theme_void(base_size = 8, base_family = "Arial") +
    theme(plot.title  = element_text(face = "bold", size = 8, hjust = 0.5),
          legend.position   = "bottom",
          legend.key.width  = unit(6, "mm"), legend.key.height = unit(1.5, "mm"),
          legend.text = element_text(size = 8),
          plot.margin = margin(2, 4, 2, 4, "pt"))
}
p5a = wrap_plots(lapply(names(a_vars), one_map), nrow = 1)
p5a

# 573. 重新拼合（B 图稍宽）并保存 ===============================
fig5 = wrap_elements(full = p5a) /
  (p5b | p5c) + plot_layout(widths = c(1.15, 1)) /
  p5d
fig5 = wrap_elements(full = p5a) / ((p5b | p5c) + plot_layout(widths = c(1.15, 1))) / p5d +
  plot_layout(heights = c(1, 1.25, 1.1)) +
  plot_annotation(tag_levels = "A")
save_jtm(fig5, "Figure5", h = 230)
# [local-only, disabled for public release] shell.exec(file.path(fd, "Figure5.pdf"))
# 573. 重新拼合（B 图稍宽）并保存 ===============================
p5bc = (p5b | p5c) + plot_layout(widths = c(1.15, 1))
fig5 = wrap_elements(full = p5a) / p5bc / p5d +
  plot_layout(heights = c(1, 1.25, 1.1)) +
  plot_annotation(tag_levels = "A")
save_jtm(fig5, "Figure5", h = 230)
# [local-only, disabled for public release] shell.exec(file.path(fd, "Figure5.pdf"))
# 574. 每个色条至少显示两个刻度 ================================
p5a = wrap_plots(lapply(names(a_vars), function(v)
  one_map(v) + scale_colour_distiller(palette = "Spectral",
                                      limits = quantile(f5$sc5[[v]], c(0.01, 0.99)),
                                      oob = scales::squish, name = NULL,
                                      breaks = scales::breaks_extended(3),
                                      labels = scales::label_number(accuracy = 0.1))),
  nrow = 1)
p5a
# 575. A 图加刻度，C 图不跟 B 图对齐，重新保存 ===================
p5a = wrap_plots(lapply(names(a_vars), function(v)
  one_map(v) + scale_colour_distiller(palette = "Spectral",
                                      limits = quantile(f5$sc5[[v]], c(0.01, 0.99)),
                                      oob = scales::squish, name = NULL,
                                      breaks = scales::breaks_extended(3),
                                      labels = scales::label_number(accuracy = 0.1))),
  nrow = 1)
p5bc = (p5b | free(p5c)) + plot_layout(widths = c(1.15, 1))
fig5 = wrap_elements(full = p5a) / p5bc / p5d +
  plot_layout(heights = c(1, 1.25, 1.1)) +
  plot_annotation(tag_levels = "A")
save_jtm(fig5, "Figure5", h = 230)
# [local-only, disabled for public release] shell.exec(file.path(fd, "Figure5.pdf"))
# 576. 改用 wrap_elements 让 C 图单独排版，再保存 ================
p5bc = (p5b | wrap_elements(full = p5c)) + plot_layout(widths = c(1.15, 1))
fig5 = wrap_elements(full = p5a) / p5bc / p5d +
  plot_layout(heights = c(1, 1.25, 1.1)) +
  plot_annotation(tag_levels = "A")
save_jtm(fig5, "Figure5", h = 230)
# [local-only, disabled for public release] shell.exec(file.path(fd, "Figure5.pdf"))
# 577. C 图 x 轴只保留 0、0.25、0.5 三个短刻度，最终保存 ==========
p5c2 = p5c + scale_x_continuous(breaks = c(0, 0.25, 0.5), labels = c("0", "0.25", "0.5"))
p5bc = (p5b | wrap_elements(full = p5c2)) + plot_layout(widths = c(1.15, 1))
fig5 = wrap_elements(full = p5a) / p5bc / p5d +
  plot_layout(heights = c(1, 1.25, 1.1)) +
  plot_annotation(tag_levels = "A")
save_jtm(fig5, "Figure5", h = 230)
# 578. 检查最终图片文件的名称和大小（每个应小于 10 MB） ==========
fd = "results/figures/final"
f  = list.files(fd, pattern = "^Figure[1-8]\\.(pdf|tiff)$", full.names = TRUE)
data.frame(file = basename(f), MB = round(file.size(f) / 1024^2, 2))
# 579. 查找补充图文件 ============================
list.files("results", pattern = "\\.(pdf|tiff|png)$",
           recursive = TRUE, full.names = TRUE)
# 582. 检查对象和列名（先运行，确认名字） ======================
library(Seurat); library(data.table)
sc_tme = readRDS("data_proc/sc_tme_scored.rds")
colnames(sc_tme@meta.data)          # 找到病人列和细胞类型列的名字
table(sc_tme$celltype)
# 如果基因集保存在文件里，在这里读入，变量名统一为 tier_sets
# tier_sets = readRDS("data_proc/tier_sets.rds")
str(tier_sets)

# 583. 设定列名（按第 582 步的结果修改） ======================
pt_col = "patient"      # 病人列
ct_col = "celltype"     # 细胞类型列
ref_ct = "Osteoblastic OS"   # 对照细胞类型（成骨型肿瘤细胞）

# 584. 做伪批量矩阵：每个“病人 × 细胞类型”取平均表达 ======
sc_tme$pb = paste(sc_tme[[pt_col]][, 1], sc_tme[[ct_col]][, 1], sep = "__")
pb = AverageExpression(sc_tme, group.by = "pb", assays = "RNA",
                       layer = "data")$RNA
pb = as.matrix(pb)
pb = pb[rowSums(pb > 0) >= 3, ]            # 去掉几乎不表达的基因
z  = t(scale(t(log1p(pb))))                # 每个基因标准化
z[is.na(z)] = 0
info = data.table(pb = colnames(z))
info[, c("pt", "ct") := tstrsplit(pb, "__")]
rm(sc_tme); gc()

# 585. 计算“某类细胞 − 肿瘤细胞”的配对差值中位数 ==============
pair_diff = function(genes, target_ct) {
  g  = intersect(genes, rownames(z))
  sc = colMeans(z[g, , drop = FALSE])
  d  = data.table(info, score = sc)
  w  = dcast(d[ct %in% c(target_ct, ref_ct)], pt ~ ct, value.var = "score")
  w  = na.omit(w)
  median(w[[target_ct]] - w[[ref_ct]])
}

# 586. 按平均表达量把基因分成 20 档，用于抽“表达量匹配”的随机基因
avg_expr = rowMeans(log1p(pb))
bins     = cut(rank(avg_expr), 20, labels = FALSE)
names(bins) = names(avg_expr)
draw_matched = function(genes) {
  g = intersect(genes, names(bins))
  unlist(lapply(g, function(x) sample(names(bins)[bins == bins[x]], 1)))
}

# 587. 对每个比较做 1000 次随机抽样 =========================
tests = data.table(
  tier = c("T1", "T2", "T3", "T4", "T4", "T4", "T1"),
  ct   = c("Osteoclast", "Myeloid", "Fibroblast/MSC",
           "Endothelial", "Pericyte", "T/NK", "Chondroblastic OS"))
set.seed(123)
res = tests[, {
  obs  = pair_diff(tier_sets[[tier]], ct)
  null = replicate(1000, pair_diff(draw_matched(tier_sets[[tier]]), ct))
  .(observed = round(obs, 3),
    null_median = round(median(null), 3),
    p_empirical = (sum(null >= obs) + 1) / 1001)
}, by = .(tier, ct)]
res
fwrite(res, "results/scRNA/random_geneset_control.csv")
