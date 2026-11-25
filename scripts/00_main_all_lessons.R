# 1. 选择磁盘：有 D 盘就用 D 盘，没有就用 C 盘
disk     = ifelse(dir.exists("D:/"), "D:/", "C:/")
root_dir = file.path(disk, "os_mn")

# 2. 定义子文件夹的名字
sub_dirs = c("data_raw", "data_proc", "scripts", "results", "ref")

# 3. 创建项目根目录
dir.create(root_dir, showWarnings = FALSE)

# 4. 创建所有子文件夹
for (d in sub_dirs) {
  dir.create(file.path(root_dir, d), showWarnings = FALSE)
}

# 5. 检查是否建好
list.dirs(root_dir, recursive = FALSE)

# 6. 把工作目录设为项目根目录
# [local-only, disabled for public release] setwd(root_dir)
getwd()
# 1. 安装读取大文件用的包
install.packages("data.table")

# 2. 安装 MR 分析的主力包（从开发者的仓库安装）
install.packages(
  "TwoSampleMR",
  repos = c("https://mrcieu.r-universe.dev", "https://cloud.r-project.org")
)

# 3. 检查是否安装成功
library(data.table)
library(TwoSampleMR)
# 1. 进入项目目录
# [local-only, disabled for public release] setwd(".")

# 2. 加载需要的包
library(data.table)

# 3. 读取锰的 GWAS 数据（文件较大，大约需要 10–60 秒）
mn_file = "data_raw/Manganese_meta.hunt_moba_pivus1.tbl"
mn      = fread(mn_file)

# 4. 看数据有多大（行数、列数）
dim(mn)

# 5. 看前 6 行
head(mn)

# 6. 看所有列名
names(mn)
# 7. 把 "P-value" 改名为 "P"（列名里有减号，在 R 里用起来不方便）
setnames(mn, "P-value", "P")

# 8. 等位基因统一转成大写（METAL 软件输出的是小写）
mn[, Allele1 := toupper(Allele1)]
mn[, Allele2 := toupper(Allele2)]

# 9. 从 MarkerName 中拆出染色体（CHR）和位置（POS）
mn[, c("CHR", "POS") := tstrsplit(sub("_.*", "", MarkerName), ":", type.convert = TRUE)]

# 10. 有多少个位点达到全基因组显著（P < 5e-8）
sum(mn$P < 5e-8)

# 11. 查找已知的锰相关位点 rs13107325（SLC39A8 基因，GRCh37 位置 4:103188709）
mn[CHR == 4 & POS == 103188709]

# 12. 保存为 R 专用格式，下次读取更快
saveRDS(mn, "data_proc/mn_gwas_clean.rds")
# ===== 第二课：筛选工具变量 =====

# 13. 进入项目目录并加载包
# [local-only, disabled for public release] setwd(".")
library(data.table)

# 14. 读取第一课保存的数据
mn = readRDS("data_proc/mn_gwas_clean.rds")

# 15. 检查有哪些染色体
table(mn$CHR)

# 16. 只保留 1–22 号常染色体，并把 CHR 转成整数
mn = mn[CHR %in% as.character(1:22)]
mn[, CHR := as.integer(CHR)]

# 17. 计算次要等位基因频率（MAF）
mn[, MAF := pmin(Freq1, 1 - Freq1)]

# 18. 按 P 值和 MAF 筛选：主分析用 5e-8，备用 5e-6
iv_main  = mn[P < 5e-8 & MAF > 0.01]
iv_relax = mn[P < 5e-6 & MAF > 0.01]
nrow(iv_main)
nrow(iv_relax)

# 19. 计算每个位点的 F 统计量
iv_main[,  F := (Effect / StdErr)^2]
iv_relax[, F := (Effect / StdErr)^2]
summary(iv_main$F)
summary(iv_relax$F)

# 20. 标记 3 个队列方向是否一致
iv_main[,  dir_ok := Direction %in% c("+++", "---")]
iv_relax[, dir_ok := Direction %in% c("+++", "---")]
table(iv_main$dir_ok)

# 21. 保存候选工具变量
fwrite(iv_main,  "data_proc/mn_iv_main.csv")
fwrite(iv_relax, "data_proc/mn_iv_relax.csv")
# ===== 第二课（补充）：设置 OpenGWAS 令牌 =====

# 22. 打开 R 的环境设置文件（.Renviron）
file.edit("~/.Renviron")
# 23. 检查令牌是否生效
library(ieugwasr)
user()
# ===== 第二课（补充）：检查令牌 =====

# 24. 检查令牌的长度
nchar(Sys.getenv("OPENGWAS_JWT"))

# 25. 检查令牌开头的 3 个字符
substr(Sys.getenv("OPENGWAS_JWT"), 1, 3)

# 26. 检查令牌里有没有空格或引号
grepl("[ \"']", Sys.getenv("OPENGWAS_JWT"))
# 28. 再次检查令牌长度（应该是几百）
nchar(Sys.getenv("OPENGWAS_JWT"))

# 29. 检查令牌是否生效
library(ieugwasr)
user()
# ===== 第二课（补充）：排查令牌问题 =====

# 30. 看 .Renviron 文件的实际位置
path.expand("~/.Renviron")

# 31. 看文件里每一行有多长（不显示内容）
nchar(readLines("~/.Renviron", warn = FALSE))

# 32. 看每一行开头的 16 个字符
substr(readLines("~/.Renviron", warn = FALSE), 1, 16)

# 33. 看 R 现在读到的令牌有多长
nchar(Sys.getenv("OPENGWAS_JWT"))
# 35. 立刻重新读取 .Renviron
readRenviron("~/.Renviron")

# 36. 再检查令牌长度
nchar(Sys.getenv("OPENGWAS_JWT"))

# 37. 检查令牌是否生效
library(ieugwasr)
user()
# ===== 第三课：转换 rsID 与 LD clumping =====

# 38. 进入项目目录并加载包
# [local-only, disabled for public release] setwd(".")
library(data.table)
library(ieugwasr)

# 39. 读取第二课保存的候选工具变量
iv = fread("data_proc/mn_iv_main.csv")

# 40. 先用已知位点测试：4:103188709 应该对应 rs13107325（SLC39A8）
variants_chrpos("4:103188709")

# 41. 生成"染色体:位置"格式的查询字符串
iv[, chrpos := paste0(CHR, ":", POS)]

# 42. 分批查询 rsID（每批 200 个，避免一次查询太大）
batches  = split(iv$chrpos, ceiling(seq_along(iv$chrpos) / 200))
res_list = lapply(batches, variants_chrpos)
rsid_tab = rbindlist(res_list, fill = TRUE)
# ===== 第三课：转换 rsID 与 LD clumping（第 43 步起，修改版）=====

# 43. 整理查询结果（这个版本没有返回基因名，只取需要的 3 列）
rsid_tab = unique(rsid_tab[, .(CHR = as.integer(chr), POS = as.integer(pos), rsid = name)])

# 44. 看看一个位置对应几个 rsID
rsid_tab[, n := .N, by = .(CHR, POS)]
table(rsid_tab$n)

# 45. 所有 rsID 都保留，合并回候选位点
iv = merge(iv, rsid_tab[, .(CHR, POS, rsid)], by = c("CHR", "POS"))
nrow(iv)

# 46. 检查：有多少个位置至少查到了一个 rsID
length(unique(iv$chrpos))

# 47. LD clumping（窗口 10000 kb，r² < 0.001，欧洲人群参考）
clump_in  = data.frame(rsid = iv$rsid, pval = iv$P, id = "manganese")
clump_out = ld_clump(clump_in, clump_kb = 10000, clump_r2 = 0.001, pop = "EUR")

# 48. 取出独立的工具变量
iv_clumped = iv[rsid %in% clump_out$rsid]
nrow(iv_clumped)

# 49. 检查：同一个位置有没有被重复保留
iv_clumped[, .N, by = .(CHR, POS)][N > 1]

# 50. 按 P 值排序，查看结果
iv_clumped[order(P), .(rsid, CHR, POS, Allele1, Allele2, Freq1, Effect, StdErr, P, F, dir_ok)]

# 51. 保存独立的工具变量
fwrite(iv_clumped, "data_proc/mn_iv_main_clumped.csv")
# 51. 保存独立的工具变量
fwrite(iv_clumped, "data_proc/mn_iv_main_clumped.csv")

# 52. 检查论文里的第 10 个位点（HSCB，22:29154455）在不在候选位点里
iv[CHR == 22 & POS == 29154455, .(rsid, P)]

# 53. 看这个位置查到了哪些 rsID
rsid_tab[CHR == 22 & POS == 29154455]

# 54. 看 rs5997397 有没有进入 clumping，以及有没有被保留下来
"rs5997397" %in% clump_in$rsid
"rs5997397" %in% clump_out$rsid
# 55. 查 HSCB 位点（3 个编号）和 TMPRSS6 位点之间的连锁程度（r）
ld = ld_matrix(c("rs5997397", "rs9620821", "rs61173979", "rs855791"), pop = "EUR", with_alleles = FALSE)
ld

# 56. 把 r 换成 r²（clumping 用的是 r²）
round(ld^2, 4)
# ===== 第四课：结局数据与孟德尔随机化 =====

# 57. 加载包和第三课的结果
# [local-only, disabled for public release] setwd(".")
library(data.table)
library(ieugwasr)
iv = fread("data_proc/mn_iv_main_clumped.csv")
nrow(iv)

# 58. 先确认令牌和额度没问题（会显示你的账号信息）
user()

# 59. 列出可能的骨肉瘤 / 骨癌数据集编号
cand = c("finn-b-C3_OSTEOSARCOMA",
         "finn-b-C3_OSTEOSARCOMA_EXALLC",
         "finn-b-C3_BONE_CARTILAGE",
         "finn-b-C3_BONE_CARTILAGE_EXALLC")

# 60. 只查这几个编号（不存在的编号会自动跳过）
os_hits = setDT(gwasinfo(cand))

# 61. 显示关键的几列
cols = intersect(c("id", "trait", "ncase", "ncontrol", "population", "year"), names(os_hits))
os_hits[, ..cols]
# 62. 安装 TwoSampleMR（只需要装一次，可能要几分钟）
install.packages("TwoSampleMR",
                 repos = c("https://mrcieu.r-universe.dev", "https://cloud.r-project.org"))
library(TwoSampleMR)

# 63. 把 9 个锰位点整理成 TwoSampleMR 要求的"暴露"格式
exp = format_data(as.data.frame(iv), type = "exposure",
                  snp_col           = "rsid",
                  beta_col          = "Effect",
                  se_col            = "StdErr",
                  effect_allele_col = "Allele1",
                  other_allele_col  = "Allele2",
                  eaf_col           = "Freq1",
                  pval_col          = "P")
exp$exposure = "Manganese"
nrow(exp)

# 64. 从 OpenGWAS 只提取这 9 个 SNP 在骨癌数据里的结果
out = extract_outcome_data(snps     = exp$SNP,
                           outcomes = "finn-b-C3_BONE_CARTILAGE_EXALLC",
                           proxies  = TRUE)
nrow(out)

# 65. 统一效应等位基因方向（harmonisation）
dat = harmonise_data(exp, out)
dat[, c("SNP", "effect_allele.exposure", "effect_allele.outcome",
        "beta.exposure", "beta.outcome", "mr_keep")]
# 66. 跑三种主要的 MR 方法
res = mr(dat, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median"))
res[, c("method", "nsnp", "b", "se", "pval")]

# 67. 把 b 换算成 OR（比值比）和 95% 置信区间
or_tab = generate_odds_ratios(res)
or_tab[, c("method", "or", "or_lci95", "or_uci95", "pval")]

# 68. 异质性检验：各个 SNP 给出的结论是否一致
mr_heterogeneity(dat)

# 69. 水平多效性检验：看 MR-Egger 的截距是否接近 0
mr_pleiotropy_test(dat)

# 70. 保存，下次可以直接读取
saveRDS(dat, "data_proc/mr_dat_bone.rds")
fwrite(or_tab, "results/mr_bone_main.csv")
# 71. 只用 2 个锰特异位点（SLC30A10 + SLC39A8）
dat_mn = subset(dat, SNP %in% c("rs6541114", "rs13107325"))
res_mn = mr(dat_mn, method_list = "mr_ivw")
generate_odds_ratios(res_mn)[, c("method", "nsnp", "or", "or_lci95", "or_uci95", "pval")]

# 72. 每个 SNP 单独算一个因果估计（Wald ratio）
ss = mr_singlesnp(dat)
ss[, c("SNP", "b", "se", "p")]

# 73. 留一法：每次去掉 1 个 SNP 重算，看结果是否被某一个位点"带偏"
loo = mr_leaveoneout(dat)
loo[, c("SNP", "b", "se", "p")]

# 74. 画散点图并保存
p1 = mr_scatter_plot(res, dat)
p1[[1]]
ggplot2::ggsave("results/mr_bone_scatter.pdf", p1[[1]], width = 6, height = 5)

# 75. 画单 SNP 森林图并保存
p2 = mr_forest_plot(ss)
p2[[1]]
ggplot2::ggsave("results/mr_bone_forest.pdf", p2[[1]], width = 6, height = 5)
loo[, c("SNP", "b", "se", "p")]
# ===== 第四课（下）：敏感性分析 =====

# 76. 安装 MR-PRESSO（只需要装一次）
install.packages("MRPRESSO",
                 repos = c("https://mrcieu.r-universe.dev", "https://cloud.r-project.org"))
library(MRPRESSO)

# 77. 跑 MR-PRESSO（随机模拟，要等一两分钟）
set.seed(123)
presso = run_mr_presso(dat, NbDistribution = 3000)
presso[[1]]$`Main MR results`
presso[[1]]$`MR-PRESSO results`$`Global Test`$Pvalue

# 78. 把 HSCB 位点（rs5997397）从第二课的结果里找回来，整理成暴露格式
hscb = fread("data_proc/mn_iv_main.csv")[CHR == 22 & POS == 29154455]
hscb[, rsid := "rs5997397"]
exp_h = format_data(as.data.frame(hscb), type = "exposure",
                    snp_col           = "rsid",
                    beta_col          = "Effect",
                    se_col            = "StdErr",
                    effect_allele_col = "Allele1",
                    other_allele_col  = "Allele2",
                    eaf_col           = "Freq1",
                    pval_col          = "P")

# 79. 提取 HSCB 在骨癌数据里的结果，统一方向，和原来 9 个合并成 10 个
out_h = extract_outcome_data(snps = "rs5997397", outcomes = "finn-b-C3_BONE_CARTILAGE_EXALLC", proxies = TRUE)
dat_h = harmonise_data(exp_h, out_h)
dat_h$id.exposure = dat$id.exposure[1]
dat_h$exposure    = "Manganese"
dat10 = as.data.frame(rbindlist(list(dat, dat_h), fill = TRUE))
nrow(dat10)

# 80. 用 10 个 SNP 重跑 MR
res10 = mr(dat10, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median"))
generate_odds_ratios(res10)[, c("method", "nsnp", "or", "or_lci95", "or_uci95", "pval")]
# ===== 第四课（下）：MVMR 调整铁 =====

# 81. 确认铁代谢的 GWAS 数据集（Benyamin 2014，4 个铁指标）
fe_info = setDT(gwasinfo(c("ieu-a-1049", "ieu-a-1050", "ieu-a-1051", "ieu-a-1052")))
fe_info[, intersect(c("id", "trait", "sample_size", "population", "year"), names(fe_info)), with = FALSE]

# 82. 提取血清铁的工具变量（P < 5e-8，已自动去除连锁）
fe_iv = extract_instruments("ieu-a-1049")
fe_iv[, c("SNP", "chr.exposure", "pos.exposure", "beta.exposure", "pval.exposure")]

# 83. 把锰和铁的 SNP 合并，每个 SNP 取两者中较小的 P 值
mv_in = rbind(data.frame(rsid = dat$SNP,   pval = dat$pval.exposure),
              data.frame(rsid = fe_iv$SNP, pval = fe_iv$pval.exposure))
mv_in = aggregate(pval ~ rsid, data = mv_in, FUN = min)
nrow(mv_in)

# 84. 对合并后的 SNP 重新做 clumping，去掉彼此连锁的
mv_in$id = "mvmr"
mv_clump = ld_clump(mv_in, clump_kb = 10000, clump_r2 = 0.001, pop = "EUR")
snps_mv = mv_clump$rsid
length(snps_mv)
snps_mv
# 85. 从完整的锰 GWAS 里取出这 10 个 SNP（读文件要等一会儿）
pos_tab = rbind(data.frame(SNP = iv$rsid,  chrpos = paste0(iv$CHR, ":", iv$POS)),
                data.frame(SNP = "rs1525892", chrpos = "3:133484712"))
mn = readRDS("data_proc/mn_gwas_clean.rds")
mn[, chrpos := paste0(CHR, ":", POS)]
mn_mv = merge(mn[chrpos %in% pos_tab$chrpos], pos_tab, by = "chrpos")
rm(mn); gc()
mn_mv[, .(SNP, chrpos, Allele1, Allele2, Effect, P)]
# 90. 看每一步还剩几个 SNP
nrow(exp_fe)
nrow(out_mv)

# 91. 分别找出：铁数据里缺的、骨癌数据里缺的、最后 MVMR 没用上的
setdiff(snps_mv, exp_fe$SNP)
setdiff(snps_mv, out_mv$SNP)
setdiff(snps_mv, rownames(mvdat$exposure_beta))

# 92. 检查回文 SNP（A/T 或 G/C），这种 SNP 在缺少频率信息时会被自动删掉
mn_mv[paste0(Allele1, Allele2) %in% c("AT", "TA", "GC", "CG"), .(SNP, Allele1, Allele2, Freq1)]

# 86. 整理成"锰"暴露格式
exp_mn = format_data(as.data.frame(mn_mv), type = "exposure",
                     snp_col = "SNP", beta_col = "Effect", se_col = "StdErr",
                     effect_allele_col = "Allele1", other_allele_col = "Allele2",
                     eaf_col = "Freq1", pval_col = "P")
exp_mn$exposure = "Manganese"; exp_mn$id.exposure = "Mn"

# 87. 从 OpenGWAS 提取这 10 个 SNP 对血清铁的效应，整理成"铁"暴露格式
fe_raw = associations(snps_mv, "ieu-a-1049", proxies = 0)
exp_fe = format_data(as.data.frame(fe_raw), type = "exposure",
                     snp_col = "rsid", beta_col = "beta", se_col = "se",
                     effect_allele_col = "ea", other_allele_col = "nea",
                     eaf_col = "eaf", pval_col = "p")
exp_fe$exposure = "Iron"; exp_fe$id.exposure = "Fe"
nrow(exp_fe)

# 88. 合并两个暴露，并提取这 10 个 SNP 的骨癌结果
mv_exp = as.data.frame(rbindlist(list(exp_mn, exp_fe), fill = TRUE))
out_mv = extract_outcome_data(snps = snps_mv, outcomes = "finn-b-C3_BONE_CARTILAGE_EXALLC", proxies = FALSE)
nrow(out_mv)

# 89. 统一方向，跑 MVMR
mvdat  = mv_harmonise_data(mv_exp, out_mv)
mv_res = mv_multiple(mvdat)
mv_res$result
# 90. 看每一步还剩几个 SNP
nrow(exp_fe)
nrow(out_mv)

# 91. 分别找出：铁数据里缺的、骨癌数据里缺的、最后 MVMR 没用上的
setdiff(snps_mv, exp_fe$SNP)
setdiff(snps_mv, out_mv$SNP)
setdiff(snps_mv, rownames(mvdat$exposure_beta))

# 92. 检查回文 SNP（A/T 或 G/C），这种 SNP 在缺少频率信息时会被自动删掉
mn_mv[paste0(Allele1, Allele2) %in% c("AT", "TA", "GC", "CG"), .(SNP, Allele1, Allele2, Freq1)]
# 93. 确认主分析里 rs6541114 是不是用代理 SNP 找到的
cols = intersect(c("SNP", "proxy.outcome", "target_snp.outcome", "proxy_snp.outcome"), names(dat))
dat[dat$SNP == "rs6541114", cols]

# 94. 重新提取铁的效应，这次允许代理 SNP，然后转换成"暴露"格式
fe_out = extract_outcome_data(snps = snps_mv, outcomes = "ieu-a-1049", proxies = TRUE)
exp_fe = convert_outcome_to_exposure(fe_out)
exp_fe$exposure = "Iron"; exp_fe$id.exposure = "Fe"
nrow(exp_fe)

# 95. 重新提取骨癌结果，同样允许代理 SNP
out_mv = extract_outcome_data(snps = snps_mv, outcomes = "finn-b-C3_BONE_CARTILAGE_EXALLC", proxies = TRUE)
nrow(out_mv)

# 96. 重新合并、统一方向、跑 MVMR
mv_exp = as.data.frame(rbindlist(list(exp_mn, exp_fe), fill = TRUE))
mvdat  = mv_harmonise_data(mv_exp, out_mv)
mv_res = mv_multiple(mvdat)
mv_res$result[, c("exposure", "nsnp", "b", "se", "pval")]
# 97. 确认铁数据里最终还缺哪个 SNP
setdiff(snps_mv, exp_fe$SNP)

# 98. 把 MVMR 结果换算成 OR 并保存
mv_tab = mv_res$result
mv_tab$or     = exp(mv_tab$b)
mv_tab$or_lci = exp(mv_tab$b - 1.96 * mv_tab$se)
mv_tab$or_uci = exp(mv_tab$b + 1.96 * mv_tab$se)
mv_tab[, c("exposure", "nsnp", "or", "or_lci", "or_uci", "pval")]
fwrite(mv_tab, "results/mvmr_bone_iron.csv")
# 99. 安装 MVMR 包（只需要装一次）
install.packages("remotes")
remotes::install_github("WSpiller/MVMR")
library(MVMR)

# 100. 把 mvdat 整理成 MVMR 包要求的格式
mv_in2 = format_mvmr(BXGs   = mvdat$exposure_beta,
                     BYG    = mvdat$outcome_beta,
                     seBXGs = mvdat$exposure_se,
                     seBYG  = mvdat$outcome_se,
                     RSID   = rownames(mvdat$exposure_beta))

# 101. 计算条件 F 统计量（每个暴露一个数字）
strength_mvmr(mv_in2, gencov = 0)

# 102. MVMR 的异质性检验（Q 统计量）
pleiotropy_mvmr(mv_in2, gencov = 0)
# ===== 第五课：单细胞分析（GSE152048） =====

# 103. 建立单细胞数据的文件夹
# [local-only, disabled for public release] setwd(".")
dir.create("data_raw/scRNA", showWarnings = FALSE)
dir.create("results/scRNA",  showWarnings = FALSE)

# 104. 查看电脑的内存有多大（单位 GB）
ram = system("powershell -command (Get-CimInstance Win32_ComputerSystem).TotalPhysicalMemory", intern = TRUE)
round(as.numeric(ram) / 1024^3, 1)

# 105. 安装 Seurat（只需要装一次，可能要 5–10 分钟）
install.packages("Seurat")

# 106. 加载 Seurat，确认版本
library(Seurat)
packageVersion("Seurat")
# 107. 准备 11 个样本的编号和下载地址
samples = c("BC2", "BC3", "BC5", "BC6", "BC10", "BC11", "BC16", "BC17", "BC20", "BC21", "BC22")
base    = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE152nnn/GSE152048/suppl/"
files   = paste0("GSE152048_", samples, ".matrix.tar.gz")
dest    = file.path("data_raw/scRNA", files)

# 108. 数据用 PowerShell 的 curl.exe --ssl-no-revoke 下载（R 的 download.file 会报 SSL 错误）

# 109. 检查下载结果（单位 MB），和 GEO 页面上的大小对比
round(file.size(dest) / 1024^2, 1)

# 110. 逐个解压到各自的文件夹
for (i in seq_along(files)) {
  untar(dest[i], exdir = file.path("data_raw/scRNA", samples[i]))
}

# 111. 看看解压出来的文件结构（以 BC2 为例）
list.files("data_raw/scRNA/BC2", recursive = TRUE)
# 112. 确认当前工作目录（应该是 .）
getwd()

# 113. 看看 scRNA 文件夹里现在有什么（有没有 .part 半截文件）
list.files("data_raw/scRNA")

# 114. 单独试下最小的一个文件（BC2，17 MB），看具体报错
download.file(paste0(base, files[1]), "data_raw/scRNA/test.tar.gz", mode = "wb")
# 115. 用 curl.exe 下载 GSE152048 的样本信息文件（很小）
mx_url  = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE152nnn/GSE152048/matrix/GSE152048_series_matrix.txt.gz"
mx_file = "data_raw/scRNA/GSE152048_series_matrix.txt.gz"
system2("curl.exe", c("-L", "--ssl-no-revoke", "-o", mx_file, mx_url))
file.size(mx_file)

# 116. 读入文件，只保留以 "!Sample_" 开头的样本信息行
mx_lines = readLines(gzfile(mx_file))
smp      = grep("^!Sample_", mx_lines, value = TRUE)

# 117. 显示样本名称和特征描述
grep("^!Sample_title|^!Sample_characteristics|^!Sample_source", smp, value = TRUE)
# 118. 挑出标题行和所有特征行
sel = grep("^!Sample_title|^!Sample_characteristics", smp, value = TRUE)

# 119. 每一行按制表符拆开、去掉引号，拼成一张表（每个 GSM 一行）
parts = lapply(sel, function(l) gsub('"', "", strsplit(l, "\t")[[1]]))
meta  = as.data.table(setNames(lapply(parts, `[`, -1), make.unique(sapply(parts, `[`, 1))))
dim(meta)

# 120. 从标题里取出病人编号（"BC2 1" → "BC2"）
meta[, patient := sub(" .*", "", `!Sample_title`)]

# 121. 每个病人只保留一行，显示所有特征
meta_pt = unique(meta[, !"!Sample_title"])
meta_pt
# 122. 根据原文给每个病人加上病灶类型
meta_pt[, lesion := "Primary"]
meta_pt[patient %in% c("BC10", "BC17"), lesion := "Lung_metastasis"]
meta_pt[patient == "BC11",              lesion := "Recurrent"]
meta_pt[patient %in% c("BC20", "BC22"), lesion := "Primary_or_Recurrent"]

# 123. 简化肿瘤类型这一列的写法
meta_pt[, histology := sub("tumor type: ", "", `!Sample_characteristics_ch1`)]
meta_pt = meta_pt[, .(patient, histology, lesion)]
meta_pt

# 124. 保存，以后读单细胞数据时用
fwrite(meta_pt, "data_proc/scRNA_sample_meta.csv")
# 109. 检查 11 个文件的大小（单位 MB）
round(file.size(dest) / 1024^2, 1)

# 110. 逐个解压到各自的文件夹
for (i in seq_along(files)) {
  untar(dest[i], exdir = file.path("data_raw/scRNA", samples[i]))
}

# 111. 看 BC2 解压出来的文件结构（显示完整路径）
list.files("data_raw/scRNA/BC2", recursive = TRUE, full.names = TRUE)
# 128. 确认 11 个样本都是同样的结构（每个应该是 3 个文件）
sapply(samples, function(s) length(list.files(file.path("data_raw/scRNA", s, s))))

# 129. 先试读 BC2，看矩阵大小（行 = 基因数，列 = 细胞数）
m = Read10X(file.path("data_raw/scRNA", "BC2", "BC2"))
dim(m)
rm(m)

# 130. 逐个读入 11 个样本，每个建成一个 Seurat 对象（要几分钟）
obj_list = lapply(samples, function(s) {
  m = Read10X(file.path("data_raw/scRNA", s, s))
  o = CreateSeuratObject(m, project = s, min.cells = 3, min.features = 200)
  o$patient = s
  o
})

# 131. 合并成一个对象，细胞编号前面加上样本名，防止重复
sc = merge(obj_list[[1]], obj_list[-1], add.cell.ids = samples)
rm(obj_list); gc()

# 132. 加上病灶类型和组织学类型，看每个样本的细胞数，保存
sc$lesion    = meta_pt$lesion[match(sc$patient, meta_pt$patient)]
sc$histology = meta_pt$histology[match(sc$patient, meta_pt$patient)]
table(sc$patient)
ncol(sc)
saveRDS(sc, "data_proc/sc_raw.rds")
# 133. 看看 BC10、BC11、BC17 多出来的是什么文件
list.files(file.path("data_raw/scRNA", c("BC10", "BC11", "BC17"), c("BC10", "BC11", "BC17")))

# 134. 确认基因名的格式，再计算每个细胞的线粒体基因比例
head(rownames(sc), 5)
sc[["percent.mt"]] = PercentageFeatureSet(sc, pattern = "^MT-")

# 135.（修正版）每个样本 3 个质控指标的中位数
qc = as.data.table(sc@meta.data)
qc[, .(cells    = .N,
       nFeature = median(as.numeric(nFeature_RNA)),
       nCount   = median(as.numeric(nCount_RNA)),
       mt_pct   = round(median(percent.mt), 1)), by = patient]

# 136. 画小提琴图，按样本分开看分布，并保存
p_qc = VlnPlot(sc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
               group.by = "patient", pt.size = 0, ncol = 1)
p_qc
ggplot2::ggsave("results/scRNA/qc_violin_before.pdf", p_qc, width = 10, height = 10)
# 137. 每个样本分别计算基因数的上限（对数尺度：中位数 + 3 倍 MAD）
qc[, f_hi := 10^(median(log10(nFeature_RNA)) + 3 * mad(log10(nFeature_RNA))), by = patient]

# 138. 标记每个细胞是否通过质控（3 个条件同时满足）
qc[, keep := nFeature_RNA >= 300 & nFeature_RNA <= f_hi & percent.mt < 20]

# 139. 看每个样本的上限，以及过滤前后的细胞数
qc[, .(f_hi   = round(f_hi[1]),
       before = .N,
       after  = sum(keep),
       kept   = paste0(round(100 * mean(keep), 1), "%")), by = patient]

# 140. 执行过滤，看剩下多少细胞
sc = sc[, qc$keep]
ncol(sc)

# 141. 保存质控后的数据（要等 1–3 分钟）
saveRDS(sc, "data_proc/sc_qc.rds")
# 142. 安装 scDblFinder（Bioconductor 的包，只需要装一次，要几分钟）
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("scDblFinder", update = FALSE, ask = FALSE)
library(scDblFinder)

# 143. 合并各样本的数据层，再转换成 scDblFinder 需要的格式
sc  = JoinLayers(sc)
sce = as.SingleCellExperiment(sc)

# 144. 按样本分别检测双细胞（可能要 10–20 分钟）
set.seed(123)
sce = scDblFinder(sce, samples = "patient")

# 145. 把结果放回 Seurat 对象，看每个样本有多少双细胞
sc$dbl = sce$scDblFinder.class
rm(sce); gc()
table(sc$patient, sc$dbl)

# 146. 只保留单细胞（singlet），保存（要等 1–3 分钟）
sc = subset(sc, dbl == "singlet")
ncol(sc)
saveRDS(sc, "data_proc/sc_qc.rds")
# 147. 标准化 → 找 2000 个高变基因 → 缩放 → PCA（要几分钟）
sc = NormalizeData(sc)
sc = FindVariableFeatures(sc, nfeatures = 2000)
sc = ScaleData(sc)
sc = RunPCA(sc, npcs = 30)

# 148. 安装并运行 Harmony，按病人校正批次效应
install.packages("harmony")
library(harmony)
sc = RunHarmony(sc, group.by.vars = "patient")

# 149. 用 Harmony 校正后的结果做聚类和 UMAP（要几分钟）
sc = FindNeighbors(sc, reduction = "harmony", dims = 1:30)
sc = FindClusters(sc, resolution = 0.5)
sc = RunUMAP(sc, reduction = "harmony", dims = 1:30)

# 150. 画两张 UMAP 图：按聚类上色、按病人上色
p1 = DimPlot(sc, group.by = "seurat_clusters", label = TRUE) + NoLegend()
p2 = DimPlot(sc, group.by = "patient")
p1 + p2
ggplot2::ggsave("results/scRNA/umap_cluster_patient.pdf", p1 + p2, width = 14, height = 6)

# 151. 保存（要等几分钟）
saveRDS(sc, "data_proc/sc_harmony.rds")
# 152. 列出各类细胞的经典标志基因
markers = list(
  Osteoblastic  = c("RUNX2", "ALPL", "IBSP", "COL1A1"),
  Chondroblastic = c("SOX9", "ACAN", "COL2A1"),
  Osteoclast    = c("ACP5", "CTSK", "MMP9"),
  Myeloid       = c("CD74", "LYZ", "CD68", "CD14"),
  T_NK          = c("CD3E", "IL7R", "NKG7", "GNLY"),
  B_Plasma      = c("MS4A1", "MZB1", "IGHG1"),
  Fibroblast    = c("DCN", "LUM"),
  MSC           = c("CXCL12", "SFRP2", "MME"),
  Endothelial   = c("PECAM1", "VWF", "EGFL7"),
  Pericyte      = c("RGS5", "ACTA2"),
  Proliferating = c("MKI67", "TOP2A")
)

# 153. 画气泡图：每个群 × 每个标志基因
p_dot = DotPlot(sc, features = markers, group.by = "seurat_clusters") + RotatedAxis()
p_dot
ggplot2::ggsave("results/scRNA/dotplot_markers.pdf", p_dot, width = 18, height = 7)

# 154. 每个群自动找出最特异的基因（每群随机抽 500 个细胞，速度快很多）
set.seed(123)
sc_ds = subset(sc, downsample = 500)
mk    = FindAllMarkers(sc_ds, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
top5  = as.data.table(mk)[order(cluster, -avg_log2FC)][, head(.SD, 5), by = cluster]
top5[, .(genes = paste(gene, collapse = ", ")), by = cluster]
rm(sc_ds); gc()
# 155. 看 17 号（肌肉）和 20 号（肺）群的细胞来自哪些病人
table(sc$seurat_clusters, sc$patient)[c("17", "20"), ]

# 156. 在 UMAP 上看浆细胞 / B 细胞标志基因的分布
p_b = FeaturePlot(sc, features = c("MZB1", "IGHG1", "MS4A1", "JCHAIN"), order = TRUE, ncol = 2)
p_b
# 156.（修正版）安装 ggrastr，再画浆细胞 / B 细胞标志基因
install.packages("ggrastr")
p_b = FeaturePlot(sc, features = c("MZB1", "MS4A1", "JCHAIN", "CD79A"), order = TRUE, ncol = 2)
p_b
# 157. 把 21 个群对应到 10 个细胞大类
anno = c("0"  = "OS_Osteoblastic",   "1"  = "Myeloid",          "2"  = "OS_Osteoblastic",
         "3"  = "OS_Osteoblastic",   "4"  = "OS_Chondroblastic", "5"  = "OS_Osteoblastic",
         "6"  = "Osteoclast",        "7"  = "T_NK",             "8"  = "Endothelial",
         "9"  = "Myeloid",           "10" = "Myeloid",          "11" = "Pericyte",
         "12" = "Myeloid",           "13" = "Fibroblast_MSC",   "14" = "OS_Osteoblastic",
         "15" = "OS_Osteoblastic",   "16" = "OS_Osteoblastic",  "17" = "Myocyte",
         "18" = "OS_Osteoblastic",   "19" = "OS_Chondroblastic", "20" = "Lung_epithelial")
sc$celltype = unname(anno[as.character(sc$seurat_clusters)])
table(sc$celltype)

# 158. 画带细胞类型标签的 UMAP
p_ct = DimPlot(sc, group.by = "celltype", label = TRUE, repel = TRUE) + NoLegend()
p_ct
ggplot2::ggsave("results/scRNA/umap_celltype.pdf", p_ct, width = 8, height = 7)

# 159. 验证：肺上皮（20 号）和肌肉（17 号）的细胞来自哪些病人
table(sc$celltype, sc$patient)[c("Lung_epithelial", "Myocyte"), ]

# 160. 保存注释好的对象（要等几分钟）
saveRDS(sc, "data_proc/sc_annotated.rds")
# 161. 定义四层锰代谢基因集
mn_sets = list(
  T1_Transport     = c("SLC39A8", "SLC39A14", "SLC30A10", "SLC11A2", "ATP2C1", "ATP13A1", "TMEM165", "TFRC"),
  T2_Enzyme        = c("SOD2", "ARG1", "ARG2", "GLUL", "PC", "PPM1A", "PPM1B"),
  T3_Glycosylation = c("B4GALT1", "GALNT1", "GALNT2", "XYLT1", "B3GAT3", "MGAT1"),
  T4_STING         = c("CGAS", "STING1", "TBK1", "IRF3", "IFNB1", "CXCL10", "ISG15")
)

# 162. 检查每个基因在数据里是否存在（不存在的会列出来）
all_g = unlist(mn_sets)
all_g[!all_g %in% rownames(sc)]

# 163. 去掉两类"污染"的正常组织细胞，只分析肿瘤及其微环境
sc_tme = subset(sc, celltype %in% c("Myocyte", "Lung_epithelial"), invert = TRUE)
ncol(sc_tme)

# 164. 画气泡图：每类细胞 × 每个锰基因
p_mn = DotPlot(sc_tme, features = mn_sets, group.by = "celltype") + RotatedAxis()
p_mn
ggplot2::ggsave("results/scRNA/dotplot_mn_genes.pdf", p_mn, width = 16, height = 6)
# 165. 确认 TMEM173 存在，看 ARG1 在原始数据里有没有
c("TMEM173", "ARG1") %in% rownames(sc)

# 166. 更新基因集：STING1 换成 TMEM173，去掉 ARG1
mn_sets$T4_STING  = c("CGAS", "TMEM173", "TBK1", "IRF3", "IFNB1", "CXCL10", "ISG15")
mn_sets$T2_Enzyme = setdiff(mn_sets$T2_Enzyme, "ARG1")

# 167. 计算两个评分：锰稳态评分（T1–T3 合并）和 STING 评分（T4）
mn_homeo = unique(unlist(mn_sets[c("T1_Transport", "T2_Enzyme", "T3_Glycosylation")]))
sc_tme   = AddModuleScore(sc_tme, features = list(mn_homeo),       name = "Mn_score")
sc_tme   = AddModuleScore(sc_tme, features = list(mn_sets$T4_STING), name = "STING_score")

# 168. 画小提琴图：两个评分在各类细胞中的分布
p_sc = VlnPlot(sc_tme, features = c("Mn_score1", "STING_score1"),
               group.by = "celltype", pt.size = 0, ncol = 1)
# 169. 四个层级分别计算评分（会生成 Tier1、Tier2、Tier3、Tier4 四列）
sc_tme = AddModuleScore(sc_tme, features = mn_sets, name = "Tier")

# 170. 每类细胞、每个评分的中位数，做成表格
sc_md = as.data.table(sc_tme@meta.data)
score_tab = sc_md[, .(T1_Transport = median(Tier1), T2_Enzyme = median(Tier2),
                      T3_Glyco     = median(Tier3), T4_STING  = median(Tier4),
                      Mn_total     = median(Mn_score1)), by = celltype]
score_tab[, (2:6) := lapply(.SD, round, 3), .SDcols = 2:6]
score_tab[order(-T1_Transport)]

# 171. 画热图：每类细胞 × 四个层级的平均评分
score_mat = as.matrix(score_tab[, .(T1_Transport, T2_Enzyme, T3_Glyco, T4_STING)])
rownames(score_mat) = score_tab$celltype
pheatmap::pheatmap(score_mat, scale = "column", cluster_cols = FALSE,
                   display_numbers = TRUE, number_format = "%.2f",
                   filename = "results/scRNA/heatmap_mn_tiers.pdf", width = 6, height = 5)
pheatmap::pheatmap(score_mat, scale = "column", cluster_cols = FALSE,
                   display_numbers = TRUE, number_format = "%.2f")
p_sc
ggplot2::ggsave("results/scRNA/vln_mn_sting_score.pdf", p_sc, width = 10, height = 8)
# 173. 关掉所有画图设备，清空状态
graphics.off()

# 174. 重新保存热图的 PDF
pheatmap::pheatmap(score_mat, scale = "column", cluster_cols = FALSE,
                   display_numbers = TRUE, number_format = "%.2f",
                   filename = "results/scRNA/heatmap_mn_tiers.pdf", width = 6, height = 5)

# 175. 检查 PDF 的大小（单位 KB），正常应该是几 KB 以上，不是 0
round(file.size("results/scRNA/heatmap_mn_tiers.pdf") / 1024, 1)

# 176. 在 Plots 窗口里再显示一次
pheatmap::pheatmap(score_mat, scale = "column", cluster_cols = FALSE,
                   display_numbers = TRUE, number_format = "%.2f")
# 172. 保存带有评分的对象（要等几分钟）
saveRDS(sc_tme, "data_proc/sc_tme_scored.rds")
fwrite(score_tab, "results/scRNA/mn_tier_scores_by_celltype.csv")
# ===== 第六课：按病人统计与病灶比较 =====

# 177. 准备：如果 R 重启过，先加载包和第五课的结果（没重启可以跳过读文件那一行）
# [local-only, disabled for public release] setwd(".")
library(data.table); library(Seurat); library(ggplot2)
# sc_tme = readRDS("data_proc/sc_tme_scored.rds")
md = as.data.table(sc_tme@meta.data)

# 178. 每个病人 × 每类细胞，计算四个层级评分的平均值（伪批量）
pb = md[, .(n_cells = .N,
            T1_Transport = mean(Tier1), T2_Enzyme = mean(Tier2),
            T3_Glyco     = mean(Tier3), T4_STING  = mean(Tier4)),
        by = .(patient, lesion, celltype)]

# 179. 只保留细胞数 ≥ 20 的组合（太少的平均值不稳定），看每类细胞有几个病人
pb = pb[n_cells >= 20]
pb[, .N, by = celltype][order(-N)]

# 180. 画箱线图：每个点是一个病人
pb_long = melt(pb, id.vars = c("patient", "lesion", "celltype", "n_cells"),
               variable.name = "tier", value.name = "score")
p_pb = ggplot(pb_long, aes(x = reorder(celltype, score), y = score)) +
  geom_boxplot(outlier.shape = NA, fill = "grey90") +
  geom_jitter(aes(color = lesion), width = 0.15, size = 1.5) +
  facet_wrap(~ tier, scales = "free_y", ncol = 2) +
  coord_flip() + theme_bw() + labs(x = NULL, y = "Mean module score per patient")
p_pb
ggsave("results/scRNA/pseudobulk_tiers_boxplot.pdf", p_pb, width = 10, height = 8)

# 181. Kruskal-Wallis 检验：每个层级的评分在不同细胞类型之间是否有差异
pb_long[, .(p_value = kruskal.test(score ~ celltype)$p.value), by = tier]
# 182. 把数据变宽：每行一个病人，每列一类细胞
wide = dcast(pb_long, patient + tier ~ celltype, value.var = "score")

# 183. 定义要比较的几组（都是和成骨型肿瘤细胞比）
comps = data.table(
  tier  = c("T1_Transport", "T1_Transport",      "T2_Enzyme", "T3_Glyco",       "T4_STING",    "T4_STING", "T4_STING"),
  cell  = c("Osteoclast",   "OS_Chondroblastic", "Myeloid",   "Fibroblast_MSC", "Endothelial", "Pericyte", "T_NK")
)

# 184. 逐组做配对 Wilcoxon 检验（只用两类细胞都有数据的病人）
comps[, c("n_pt", "median_diff", "p") := {
  d  = wide[tier == .BY$tier]
  x  = d[[cell]]; y = d[["OS_Osteoblastic"]]
  ok = !is.na(x) & !is.na(y)
  list(sum(ok), round(median(x[ok] - y[ok]), 3), wilcox.test(x[ok], y[ok], paired = TRUE)$p.value)
}, by = .(tier, cell)]

# 185. 多重检验校正（BH 法），显示结果
comps[, p_adj := p.adjust(p, method = "BH")]
comps[, `:=`(p = signif(p, 3), p_adj = signif(p_adj, 3))]
comps
fwrite(comps, "results/scRNA/pseudobulk_paired_tests.csv")
# 185b. 统计每组比较中，有几个病人的差值是正的（方向一致）
comps[, n_consistent := {
  d = wide[tier == .BY$tier]; x = d[[cell]]; y = d[["OS_Osteoblastic"]]
  ok = !is.na(x) & !is.na(y); sum(x[ok] > y[ok])
}, by = .(tier, cell)]
comps[, .(tier, cell, n_pt, n_consistent, p_adj)]
fwrite(comps, "results/scRNA/pseudobulk_paired_tests.csv")
# 186. 只看成骨型肿瘤细胞，把病灶分成两组：原发 vs 转移/复发（去掉分类不明的 BC20、BC22）
os_pb = pb[celltype == "OS_Osteoblastic" & lesion != "Primary_or_Recurrent"]
os_pb[, group := ifelse(lesion == "Primary", "Primary", "Met_or_Rec")]
os_pb[, .(patient, lesion, group, T1_Transport, T2_Enzyme, T3_Glyco, T4_STING)]

# 187. 两组的中位数和 Wilcoxon 检验（样本量小，仅供参考）
os_long = melt(os_pb, id.vars = c("patient", "group"),
               measure.vars = c("T1_Transport", "T2_Enzyme", "T3_Glyco", "T4_STING"),
               variable.name = "tier", value.name = "score")
os_long[, .(Primary    = round(median(score[group == "Primary"]), 3),
            Met_or_Rec = round(median(score[group == "Met_or_Rec"]), 3),
            p          = signif(wilcox.test(score ~ group)$p.value, 3)), by = tier]

# 188. 画图：每个点是一个病人
p_les = ggplot(os_long, aes(x = group, y = score, color = group)) +
  geom_boxplot(outlier.shape = NA) + geom_point(size = 2) +
  geom_text(aes(label = patient), size = 3, vjust = -0.8, color = "black") +
  facet_wrap(~ tier, scales = "free_y", ncol = 4) + theme_bw() + NoLegend() +
  labs(x = NULL, y = "Mean module score (OS osteoblastic cells)")
p_les
ggsave("results/scRNA/os_lesion_compare.pdf", p_les, width = 12, height = 4)
# ===== 第七课：细胞通讯（CellChat） =====

# 189. 安装 CellChat 及其依赖（只需要装一次，可能要 10 分钟以上）
BiocManager::install(c("ComplexHeatmap", "BiocNeighbors"), update = FALSE, ask = FALSE)
install.packages("CellChat", repos = c("https://jinworks.r-universe.dev", "https://cloud.r-project.org"))
library(CellChat)
packageVersion("CellChat")

# 190. 每类细胞最多随机抽 2000 个
Idents(sc_tme) = "celltype"
set.seed(123)
sc_cc = subset(sc_tme, downsample = 2000)
table(sc_cc$celltype)

# 191. 建立 CellChat 对象，使用人类配体–受体数据库
cellchat = createCellChat(object = sc_cc, group.by = "celltype", assay = "RNA")
cellchat@DB = subsetDB(CellChatDB.human)
rm(sc_cc); gc()

# 192. 预处理：找出每类细胞高表达的配体和受体
cellchat = subsetData(cellchat)
cellchat = identifyOverExpressedGenes(cellchat)
cellchat = identifyOverExpressedInteractions(cellchat)
# 189.（修正版）从 blaserlab 的 r-universe 安装 CellChat（Windows 编译好的版本）
install.packages("CellChat", repos = c("https://blaserlab.r-universe.dev", "https://cloud.r-project.org"))
library(CellChat)
packageVersion("CellChat")
# 193. 检查几个关键依赖是否已经装好
c("NMF", "ComplexHeatmap", "BiocNeighbors", "circlize", "ggalluvial") %in% rownames(installed.packages())

file.exists("data_proc/sc_tme_scored.rds")
saveRDS(sc_tme, "data_proc/sc_tme_scored.rds")
rm(list = ls()); gc()

# 198. 用 curl.exe 下载 CellChat 的 Windows 安装包
cc_url  = "https://blaserlab.r-universe.dev/bin/windows/contrib/4.6/CellChat_2.2.0.9001.zip"
cc_file = "ref/CellChat_2.2.0.9001.zip"
system2("curl.exe", c("-L", "--ssl-no-revoke", "-o", cc_file, cc_url))
file.size(cc_file)
# 200. 改用清华镜像，并确认 Rtools 能用（应该显示 TRUE）
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
install.packages("pkgbuild")
pkgbuild::has_build_tools(debug = TRUE)

# 201. 从 GitHub 下载 CellChat 源代码并编译安装（依赖包会从清华镜像下载）
install.packages("remotes")
remotes::install_github("jinworks/CellChat", upgrade = "never")
library(CellChat)
packageVersion("CellChat")
# 197. 重新加载包和第六课保存的数据（读文件要等一会儿）
# [local-only, disabled for public release] setwd(".")
library(data.table); library(Seurat); library(ggplot2); library(CellChat)
sc_tme = readRDS("data_proc/sc_tme_scored.rds")
ncol(sc_tme)

# 190. 每类细胞最多随机抽 2000 个
Idents(sc_tme) = "celltype"
set.seed(123)
sc_cc = subset(sc_tme, downsample = 2000)
table(sc_cc$celltype)

# 191. 建立 CellChat 对象，使用人类配体–受体数据库
cellchat = createCellChat(object = sc_cc, group.by = "celltype", assay = "RNA")
cellchat@DB = subsetDB(CellChatDB.human)
rm(sc_cc); gc()

# 192. 预处理：找出每类细胞高表达的配体和受体
cellchat = subsetData(cellchat)
cellchat = identifyOverExpressedGenes(cellchat)
cellchat = identifyOverExpressedInteractions(cellchat)
# 192b. 安装 presto（加速 Wilcoxon 检验，只需要装一次）
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
remotes::install_github("immunogenomics/presto", upgrade = "never")

# 192.（续）重新运行报错的那两行
cellchat = identifyOverExpressedGenes(cellchat)
cellchat = identifyOverExpressedInteractions(cellchat)
# 202. 计算每对细胞之间、每个配体–受体的通讯概率（最慢的一步，可能要 10–30 分钟）
cellchat = computeCommunProb(cellchat)

# 203. 过滤 + 汇总：去掉细胞太少的通讯，按信号通路汇总，再算整体网络
cellchat = filterCommunication(cellchat, min.cells = 10)
cellchat = computeCommunProbPathway(cellchat)
cellchat = aggregateNet(cellchat)

# 204. 保存 CellChat 对象
saveRDS(cellchat, "data_proc/cellchat.rds")

# 205. 画两张圆圈图：细胞之间通讯的"次数"和"强度"，并保存
groupSize = as.numeric(table(cellchat@idents))
pdf("results/scRNA/cellchat_circle.pdf", width = 12, height = 6)
par(mfrow = c(1, 2), xpd = TRUE)
netVisual_circle(cellchat@net$count,  vertex.weight = groupSize, weight.scale = TRUE,
                 label.edge = FALSE, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = TRUE,
                 label.edge = FALSE, title.name = "Interaction strength")
dev.off()
# 206. 列出所有显著的信号通路
cellchat@netP$pathways

# 207. 提取"锰活跃细胞 → 肿瘤细胞"的所有配体–受体对，按通讯概率排序
mn_src = c("Osteoclast", "Myeloid", "Fibroblast_MSC")
os_tgt = c("OS_Osteoblastic", "OS_Chondroblastic")
lr = as.data.table(subsetCommunication(cellchat, sources.use = mn_src, targets.use = os_tgt))
lr[order(-prob)][1:20, .(source, target, ligand, receptor, pathway_name, prob = signif(prob, 3))]
fwrite(lr, "results/scRNA/cellchat_mn_to_os.csv")

# 208. 每类锰活跃细胞贡献了哪些主要通路（按总概率排）
lr[, .(n_pairs = .N, total_prob = signif(sum(prob), 3)), by = .(source, pathway_name)][order(source, -total_prob)][, head(.SD, 5), by = source]

# 209. 画气泡图：锰活跃细胞 → 肿瘤细胞的主要配体–受体对
p_bub = netVisual_bubble(cellchat, sources.use = mn_src, targets.use = os_tgt,
                         remove.isolate = TRUE, angle.x = 45)
ggsave("results/scRNA/cellchat_bubble_mn_to_os.pdf", p_bub, width = 10, height = 14)
# 210. 提取全部细胞间通讯，挑出受体是整合素（名字里带 ITG）的
all_lr = as.data.table(subsetCommunication(cellchat))
itg = all_lr[grepl("ITG", receptor)]
nrow(itg)

# 211. 整合素信号按"发出者 × 接收者"汇总，画热图
itg_mat = dcast(itg[, .(prob = sum(prob)), by = .(source, target)],
                source ~ target, value.var = "prob", fill = 0)
itg_m   = as.matrix(itg_mat[, -1]); rownames(itg_m) = itg_mat$source
pheatmap::pheatmap(itg_m, cluster_rows = FALSE, cluster_cols = FALSE,
                   display_numbers = TRUE, number_format = "%.2f",
                   main = "Integrin signalling (sum of probability)",
                   filename = "results/scRNA/cellchat_integrin_heatmap.pdf", width = 7, height = 6)
pheatmap::pheatmap(itg_m, cluster_rows = FALSE, cluster_cols = FALSE,
                   display_numbers = TRUE, number_format = "%.2f",
                   main = "Integrin signalling (sum of probability)")

# 212. 看关键受体在各类细胞上的表达（整合素、CD44、BSG、PLXNB1）
rec = c("ITGA10", "ITGA6", "ITGB1", "ITGAV", "ITGB3", "CD44", "BSG", "PLXNB1", "NCL", "SDC2")
p_rec = DotPlot(sc_tme, features = rec, group.by = "celltype") + RotatedAxis()
p_rec
ggsave("results/scRNA/dotplot_receptors.pdf", p_rec, width = 9, height = 5)

# 213. 反向：肿瘤细胞 → 破骨细胞 / 髓系细胞，通讯概率最高的 15 对
rev_lr = all_lr[source %in% os_tgt & target %in% c("Osteoclast", "Myeloid")]
rev_lr[order(-prob)][1:15, .(source, target, ligand, receptor, pathway_name, prob = signif(prob, 3))]
fwrite(rev_lr, "results/scRNA/cellchat_os_to_mn.csv")
# ===== 第八课：TARGET-OS 队列与生存分析 =====

# 214. 建文件夹，用 curl.exe 下载 4 个文件（表达、基因对照表、临床、生存）
dir.create("data_raw/bulk", showWarnings = FALSE)
xena  = "https://gdc.xenahubs.net/download/"
files = c("TARGET-OS.star_tpm.tsv.gz", "gencode.v36.annotation.gtf.gene.probemap",
          "TARGET-OS.clinical.tsv.gz",  "TARGET-OS.survival.tsv.gz")
for (f in files) {
  system2("curl.exe", c("-L", "--ssl-no-revoke", "-o", file.path("data_raw/bulk", f), paste0(xena, f)))
}

# 215. 检查文件大小（单位 MB）
round(file.size(file.path("data_raw/bulk", files)) / 1024^2, 2)

# 216. 读入表达矩阵和基因对照表，看看前几行
expr  = fread("data_raw/bulk/TARGET-OS.star_tpm.tsv.gz")
pmap  = fread("data_raw/bulk/gencode.v36.annotation.gtf.gene.probemap")
dim(expr)
expr[1:3, 1:4]
head(pmap, 3)

# 217. 读入临床和生存数据，看有哪些列
clin = fread("data_raw/bulk/TARGET-OS.clinical.tsv.gz")
surv = fread("data_raw/bulk/TARGET-OS.survival.tsv.gz")
dim(clin); dim(surv)
names(clin)
head(surv, 3)
# 218. 表达矩阵：Ensembl 编号换成基因名，同名基因取平均，变成"基因 × 样本"的矩阵
expr = merge(pmap[, .(Ensembl_ID = id, gene)], expr, by = "Ensembl_ID")
expr[, Ensembl_ID := NULL]
expr = expr[, lapply(.SD, mean), by = gene]
mat  = as.matrix(expr[, -1]); rownames(mat) = expr$gene
dim(mat)

# 219. 看样本类型：条形码第 18–19 位，01 = 原发肿瘤，其他是复发/转移等
table(substr(colnames(mat), 18, 19))

# 220. 整理临床表：每个样本一行，自己计算生存时间和状态
cl = clin[sample %in% colnames(mat),
          .(sample,
            status = vital_status.demographic,
            d_death = days_to_death.demographic,
            d_fu    = days_to_last_follow_up.diagnoses,
            meta    = metastasis_at_diagnosis.diagnoses,
            age     = age_at_diagnosis.diagnoses / 365.25,
            sex     = gender.demographic)]
cl = unique(cl, by = "sample")
cl[, OS.time := ifelse(status == "Dead", d_death, d_fu)]
cl[, OS      := as.integer(status == "Dead")]
nrow(cl)

# 221. 看看各项信息的分布
table(cl$status, useNA = "ifany")
table(cl$meta,   useNA = "ifany")
summary(cl$OS.time)
sum(!is.na(cl$OS.time) & !is.na(cl$OS))
# 222. 看一个样本有几行，随访时间在不同行里是否不一样
clin[sample %in% colnames(mat), .N, by = sample][, table(N)]
clin[sample %in% colnames(mat), .(n_fu = sum(!is.na(days_to_last_follow_up.diagnoses))), by = sample][, table(n_fu)]

# 223. 每个样本把所有行合起来：取非空的最大值（随访时间、死亡时间）
max_na = function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
fu = clin[sample %in% colnames(mat),
          .(d_fu2    = max_na(days_to_last_follow_up.diagnoses),
            d_death2 = max_na(days_to_death.demographic)), by = sample]
cl = merge(cl, fu, by = "sample")
cl[, OS.time := ifelse(status == "Dead", d_death2, d_fu2)]

# 224. 再用 Xena 的生存文件补充：自己算的是 NA 时，用 Xena 的
cl = merge(cl, surv[, .(sample, xena_time = OS.time, xena_os = OS)], by = "sample", all.x = TRUE)
cl[is.na(OS.time) & !is.na(xena_time), `:=`(OS.time = xena_time, OS = xena_os)]
cl[status == "Not Reported" & is.na(xena_os), OS := NA]

# 225. 最终有多少人生存数据完整，其中多少人去世
cl[!is.na(OS.time) & !is.na(OS), .(n = .N, dead = sum(OS), median_days = median(OS.time))]
# 226. 用 GDC 接口查询 TARGET-OS 的"临床补充文件"列表
filt = '{"op":"and","content":[{"op":"in","content":{"field":"cases.project.project_id","value":["TARGET-OS"]}},{"op":"in","content":{"field":"data_type","value":["Clinical Supplement"]}}]}'
q_url = paste0("https://api.gdc.cancer.gov/files?size=50&fields=file_id,file_name,file_size,access&filters=",
               URLencode(filt, reserved = TRUE))
system2("curl.exe", c("-s", "--ssl-no-revoke", "-o", "data_raw/bulk/gdc_clin_files.json", shQuote(q_url, type = "cmd")))

# 227. 读取查询结果，显示文件列表
gdc_files = jsonlite::fromJSON("data_raw/bulk/gdc_clin_files.json")$data$hits
gdc_files
# 228. 下载第 3、4 个文件（发现队列和验证队列的临床数据）
dl = function(id, out) system2("curl.exe", c("-s", "--ssl-no-revoke", "-o", out,
                                             paste0("https://api.gdc.cancer.gov/data/", id)))
dl("0bbb355b-da04-474a-a4f8-cacac10af5ee", "data_raw/bulk/TARGET_OS_ClinicalData_Discovery.xlsx")
dl("7a862330-297a-466c-bca5-6c4318c286e6", "data_raw/bulk/TARGET_OS_ClinicalData_Validation.xlsx")
file.size(c("data_raw/bulk/TARGET_OS_ClinicalData_Discovery.xlsx",
            "data_raw/bulk/TARGET_OS_ClinicalData_Validation.xlsx"))

# 229. 读入两个 Excel 文件，合并，看列名
library(readxl)
cd = rbind(as.data.table(read_excel("data_raw/bulk/TARGET_OS_ClinicalData_Discovery.xlsx")),
           as.data.table(read_excel("data_raw/bulk/TARGET_OS_ClinicalData_Validation.xlsx")),
           fill = TRUE)
dim(cd)
names(cd)

# 230. 看前 3 行
head(cd, 3)
# 231. 从原始临床表里挑出需要的列，改成简短的名字（同一个病人有多行时保留第一行）
cd2 = cd[, .(USI         = `TARGET USI`,
             vital       = `Vital Status`,
             os_days     = as.numeric(`Overall Survival Time in Days`),
             efs_days    = as.numeric(`Event Free Survival Time in Days`),
             first_event = `First Event`,
             dz_dx       = `Disease at diagnosis`,
             age_y       = `Age at Diagnosis in Days` / 365.25,
             sex         = Gender,
             hist_resp   = `Histologic response`)]
cd2 = unique(cd2, by = "USI")

# 232. 表达样本名的前 16 位就是病人编号，和临床表合并
cl2 = data.table(sample = colnames(mat), USI = substr(colnames(mat), 1, 16))
cl2 = merge(cl2, cd2, by = "USI", all.x = TRUE)
nrow(cl2)

# 233. 看看几个关键变量的取值
table(cl2$vital,       useNA = "ifany")
table(cl2$dz_dx,       useNA = "ifany")
table(cl2$first_event, useNA = "ifany")

# 234. 总生存数据完整的人数、去世人数、随访中位数
cl2[, OS := as.integer(vital == "Dead")]
cl2[!is.na(os_days) & vital %in% c("Alive", "Dead"),
    .(n = .N, dead = sum(OS), median_fu_days = median(os_days))]
# 235. 定义四层基因集（bulk 用新版名称 STING1），检查哪些基因不在数据里
bulk_sets = list(
  T1_Transport     = c("SLC39A8", "SLC39A14", "SLC30A10", "SLC11A2", "ATP2C1", "ATP13A1", "TMEM165", "TFRC"),
  T2_Enzyme        = c("SOD2", "ARG2", "GLUL", "PC", "PPM1A", "PPM1B"),
  T3_Glycosylation = c("B4GALT1", "GALNT1", "GALNT2", "XYLT1", "B3GAT3", "MGAT1"),
  T4_STING         = c("CGAS", "STING1", "TBK1", "IRF3", "IFNB1", "CXCL10", "ISG15")
)
setdiff(unlist(bulk_sets), rownames(mat))

# 236. 每个基因在样本间标准化（z 分数），再按层级取平均，得到每个样本 4 个评分
z = t(scale(t(mat[unlist(bulk_sets), ])))
bs = as.data.table(sapply(bulk_sets, function(g) colMeans(z[g, , drop = FALSE])), keep.rownames = "sample")
head(bs, 3)

# 237. 和临床表合并，整理转移分组和无事件生存
bk = merge(cl2, bs, by = "sample")
bk[, meta := fifelse(grepl("^Metastatic", dz_dx), "Metastatic",
                     fifelse(grepl("^Non", dz_dx), "Non-metastatic", NA_character_))]
bk[, EFS := as.integer(first_event %in% c("Relapse", "Death", "Second Malignant Neoplasm"))]
table(bk$meta, useNA = "ifany")

# 238. 转移 vs 未转移：四个评分的比较（Wilcoxon 检验）
bk_long = melt(bk[!is.na(meta)], id.vars = c("sample", "meta"),
               measure.vars = names(bulk_sets), variable.name = "tier", value.name = "score")
bk_long[, .(Non_met = round(median(score[meta == "Non-metastatic"]), 3),
            Met     = round(median(score[meta == "Metastatic"]), 3),
            p       = signif(wilcox.test(score ~ meta)$p.value, 3)), by = tier]

# 239. 单因素 Cox 回归：每个评分（标准化后）和总生存、无事件生存的关系
library(survival)
cox_one = function(tier, time, event) {
  d = bk[!is.na(get(time)) & !is.na(get(event))]
  f = as.formula(paste0("Surv(", time, ", ", event, ") ~ scale(", tier, ")"))
  s = summary(coxph(f, data = d))
  data.table(tier, endpoint = event, n = nrow(d),
             HR = round(s$conf.int[1], 2), lo = round(s$conf.int[3], 2),
             hi = round(s$conf.int[4], 2), p = signif(s$coefficients[5], 3))
}
rbindlist(c(lapply(names(bulk_sets), cox_one, time = "os_days",  event = "OS"),
            lapply(names(bulk_sets), cox_one, time = "efs_days", event = "EFS")))
# 240. bulk 的锰评分是否主要反映细胞组成：和各类细胞标志基因的相关性
mk_bulk = c(Osteoclast = "ACP5", Osteoclast2 = "CTSK", Myeloid = "CD68", Fibroblast = "DCN", T_cell = "CD3E", Endothelial = "PECAM1")
cor_tab = sapply(names(bulk_sets), function(t) round(sapply(mk_bulk, function(g) cor(bs[[t]], mat[g, bs$sample], method = "spearman")), 2))
cor_tab
# 241. 保存 TARGET 的表达矩阵、临床表和评分
saveRDS(mat, "data_proc/target_expr.rds")
fwrite(bk,   "data_proc/target_clin_scores.csv")
# ===== 第九课：外部验证（GSE21257） =====

# 242. 下载表达矩阵（series matrix）和芯片注释文件（探针 → 基因名）
dir.create("data_raw/gse21257", showWarnings = FALSE)
g_mx  = "data_raw/gse21257/GSE21257_series_matrix.txt.gz"
g_ann = "data_raw/gse21257/GPL6102.annot.gz"
system2("curl.exe", c("-L", "--ssl-no-revoke", "-o", g_mx,
                      "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE21nnn/GSE21257/matrix/GSE21257_series_matrix.txt.gz"))
system2("curl.exe", c("-L", "--ssl-no-revoke", "-o", g_ann,
                      "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL6nnn/GPL6102/annot/GPL6102.annot.gz"))
round(file.size(c(g_mx, g_ann)) / 1024^2, 2)

# 243. 读入 series matrix，整理样本的临床信息（每个样本一行）
mx    = readLines(gzfile(g_mx))
ch    = grep("^!Sample_geo_accession|^!Sample_characteristics_ch1", mx, value = TRUE)
parts = lapply(ch, function(l) gsub('"', "", strsplit(l, "\t")[[1]][-1]))
meta21 = data.table(sample = parts[[1]])
for (v in parts[-1]) meta21[[sub(":.*", "", v[1])]] = sub("^[^:]*: ", "", v)
dim(meta21)
names(meta21)

# 244. 看每一列有哪些取值（数字很多的列只显示前几个）
lapply(meta21[, -1], function(x) head(sort(table(x), decreasing = TRUE), 6))

# 245. 读入表达矩阵部分，看大小和数值范围
b  = grep("^!series_matrix_table_begin", mx)
e  = grep("^!series_matrix_table_end",   mx)
ex = fread(text = mx[(b + 1):(e - 1)])
dim(ex)
summary(as.numeric(as.matrix(ex[, -1])))
# 246. 从文字里拆出生存状态、生存时间（月）和转移信息
meta21[, OS      := as.integer(grepl("^Deceased", status))]
meta21[, os_mon  := as.numeric(sub(".* at ([0-9]+) months.*", "\\1", status))]
meta21[, meta_any := fifelse(group == "No metastases", "No", "Yes")]
meta21[, meta_dx  := fifelse(group == "Metastases present at diagnosis", "Yes", "No")]
meta21[, .(n = .N, dead = sum(OS), median_mon = median(os_mon)), by = meta_any]

# 247. 读入芯片注释，探针换成基因名；同一基因有多个探针时，保留平均表达最高的那个
an  = readLines(gzfile(g_ann))
ab  = grep("^!platform_table_begin", an); ae = grep("^!platform_table_end", an)
ann = fread(text = an[(ab + 1):(ae - 1)], select = c("ID", "Gene symbol"))
setnames(ann, c("ID_REF", "gene"))
ex2 = merge(ann[gene != ""], ex, by = "ID_REF")
ex2[, avg := rowMeans(.SD), .SDcols = meta21$sample]
ex2 = ex2[order(-avg)][!duplicated(gene)]
mat21 = as.matrix(ex2[, meta21$sample, with = FALSE]); rownames(mat21) = ex2$gene
dim(mat21)
setdiff(unlist(bulk_sets), rownames(mat21))

# 248. 用和 TARGET 完全一样的方法算四个评分（z 分数平均），只用芯片上有的基因
sets21 = lapply(bulk_sets, intersect, rownames(mat21))
z21    = t(scale(t(mat21[unlist(sets21), ])))
bs21   = as.data.table(sapply(sets21, function(g) colMeans(z21[g, , drop = FALSE])), keep.rownames = "sample")
g21    = merge(meta21, bs21, by = "sample")

# 249. 转移（有 vs 无）的比较
g21_long = melt(g21, id.vars = c("sample", "meta_any"), measure.vars = names(bulk_sets),
                variable.name = "tier", value.name = "score")
g21_long[, .(No  = round(median(score[meta_any == "No"]),  3),
             Yes = round(median(score[meta_any == "Yes"]), 3),
             p   = signif(wilcox.test(score ~ meta_any)$p.value, 3)), by = tier]

# 250. 单因素 Cox 回归：每个评分和总生存的关系
rbindlist(lapply(names(bulk_sets), function(t) {
  s = summary(coxph(as.formula(paste0("Surv(os_mon, OS) ~ scale(", t, ")")), data = g21))
  data.table(tier = t, HR = round(s$conf.int[1], 2), lo = round(s$conf.int[3], 2),
             hi = round(s$conf.int[4], 2), p = signif(s$coefficients[5], 3))
}))
# 249a. 确认问题：看 meta_any 的分布，以及 group 原文前后有没有空格
table(g21$meta_any)
unique(sprintf("[%s]", meta21$group))

# 246.（修正版）用 grepl 匹配关键词，不受空格影响；然后重新合并
meta21[, meta_any := fifelse(grepl("No metastases", group), "No", "Yes")]
meta21[, meta_dx  := fifelse(grepl("present at diagnosis", group), "Yes", "No")]
table(meta21$meta_any)
g21 = merge(meta21, bs21, by = "sample")
# 249.（修正版）转移（有 vs 无）的比较：不用公式，直接传两组数值
g21_long = melt(g21, id.vars = c("sample", "meta_any"), measure.vars = names(bulk_sets),
                variable.name = "tier", value.name = "score")
g21_long[, .(No  = round(median(score[meta_any == "No"]),  3),
             Yes = round(median(score[meta_any == "Yes"]), 3),
             p   = signif(wilcox.test(score[meta_any == "Yes"], score[meta_any == "No"])$p.value, 3)),
         by = tier]

# 251. 诊断：每个层级在芯片上有几个基因；哪些基因的 z 分数是 NaN；哪个评分是 NA
lengths(sets21)
names(which(rowSums(is.na(z21)) > 0))
colSums(is.na(bs21[, -1]))

# 252.（修正）去掉在样本间没有变化的基因，重新计算评分
keep_g = rownames(z21)[rowSums(is.na(z21)) == 0]
sets21 = lapply(sets21, intersect, keep_g)
lengths(sets21)
bs21   = as.data.table(sapply(sets21, function(g) colMeans(z21[g, , drop = FALSE])), keep.rownames = "sample")
g21    = merge(meta21, bs21, by = "sample")
colSums(is.na(bs21[, -1]))
# 253. 诊断第 247 步：注释表、表达矩阵的探针编号、合并后的结果
dim(ann); head(ann, 3)
head(ex$ID_REF, 3)
sum(ex$ID_REF %in% ann$ID_REF)
dim(ex2)
dim(mat21); head(rownames(mat21), 5)

# 254. 安装 nuID 注释包（Bioconductor，只需要装一次）
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
BiocManager::install("lumiHumanAll.db", update = FALSE, ask = FALSE)
library(lumiHumanAll.db)

# 255. 用 nuID 查基因名，看能匹配上多少
sym = AnnotationDbi::mapIds(lumiHumanAll.db, keys = ex$ID_REF, column = "SYMBOL",
                            keytype = "PROBEID", multiVals = "first")
sum(!is.na(sym))

# 256. 重新做第 247 步：nuID 换成基因名，同一基因保留平均表达最高的探针
ann = data.table(ID_REF = names(sym), gene = unname(sym))[!is.na(gene)]
ex2 = merge(ann, ex, by = "ID_REF")
ex2[, avg := rowMeans(.SD), .SDcols = meta21$sample]
ex2 = ex2[order(-avg)][!duplicated(gene)]
mat21 = as.matrix(ex2[, meta21$sample, with = FALSE]); rownames(mat21) = ex2$gene
dim(mat21)
setdiff(unlist(bulk_sets), rownames(mat21))
# 248. 用和 TARGET 完全一样的方法算四个评分
sets21 = lapply(bulk_sets, intersect, rownames(mat21))
z21    = t(scale(t(mat21[unlist(sets21), ])))
bs21   = as.data.table(sapply(sets21, function(g) colMeans(z21[g, , drop = FALSE])), keep.rownames = "sample")
g21    = merge(meta21, bs21, by = "sample")
colSums(is.na(bs21[, -1]))

# 249. 转移（有 vs 无）的比较
g21_long = melt(g21, id.vars = c("sample", "meta_any"), measure.vars = names(bulk_sets),
                variable.name = "tier", value.name = "score")
g21_long[, .(No  = round(median(score[meta_any == "No"]),  3),
             Yes = round(median(score[meta_any == "Yes"]), 3),
             p   = signif(wilcox.test(score[meta_any == "Yes"], score[meta_any == "No"])$p.value, 3)),
         by = tier]

# 250. 单因素 Cox 回归：每个评分和总生存的关系
rbindlist(lapply(names(bulk_sets), function(t) {
  s = summary(coxph(as.formula(paste0("Surv(os_mon, OS) ~ scale(", t, ")")), data = g21))
  data.table(tier = t, HR = round(s$conf.int[1], 2), lo = round(s$conf.int[3], 2),
             hi = round(s$conf.int[4], 2), p = signif(s$coefficients[5], 3))
}))
# 257. GSE21257：STING 评分是否也主要反映 T 细胞、内皮细胞浸润（和 TARGET 对照）
round(sapply(c(T_cell = "CD3E", Endothelial = "PECAM1", Myeloid = "CD68", Fibroblast = "DCN"),
             function(g) cor(g21$T4_STING, mat21[g, g21$sample], method = "spearman")), 2)

# 258. 多因素 Cox：T3 糖基化在调整 Huvos 分级、诊断时转移、年龄后是否仍有意义
g21[, huvos_good := fifelse(`huvos grade` %in% c("3", "4"), 1L, fifelse(`huvos grade` %in% c("1", "2"), 0L, NA_integer_))]
g21[, age_y := as.numeric(sub(" months", "", age)) / 12]
summary(coxph(Surv(os_mon, OS) ~ scale(T3_Glycosylation) + huvos_good + meta_dx + age_y, data = g21))$conf.int

# 259. 保存 GSE21257 的结果
fwrite(g21, "data_proc/gse21257_clin_scores.csv")
# ===== 第十课：随机基因集对照与免疫浸润 =====

# 260. GSE21257：所有基因都做 z 分数标准化（去掉有缺失的），准备好转移和生存信息
z_all = t(scale(t(mat21)))
z_all = z_all[complete.cases(z_all), ]
idx   = match(colnames(z_all), g21$sample)
yes   = g21$meta_any[idx] == "Yes"
tm    = g21$os_mon[idx]; ev = g21$OS[idx]
dim(z_all)

# 261. 两个检验函数：转移比较的 P 值、生存 Cox 的 P 值
p_meta = function(s) wilcox.test(s[yes], s[!yes])$p.value
p_cox  = function(s) summary(coxph(Surv(tm, ev) ~ scale(s)))$coefficients[5]

# 262. 每个层级：随机抽同样数量的基因 1000 次，算经验 P 值（要 1–3 分钟）
set.seed(123)
B = 1000
emp = rbindlist(lapply(names(sets21), function(t) {
  n   = length(sets21[[t]])
  obs = colMeans(z_all[sets21[[t]], ])
  rnd = replicate(B, colMeans(z_all[sample(nrow(z_all), n), ]))
  rp_meta = apply(rnd, 2, p_meta)
  rp_cox  = apply(rnd, 2, p_cox)
  data.table(tier = t, n_genes = n,
             p_meta = signif(p_meta(obs), 3), emp_meta = mean(rp_meta <= p_meta(obs)),
             p_cox  = signif(p_cox(obs), 3),  emp_cox  = mean(rp_cox  <= p_cox(obs)))
}))
emp
fwrite(emp, "results/random_geneset_gse21257.csv")
# 263. 用经典标志基因算几类免疫/间质细胞的浸润评分（不含 STING 基因，避免重复）
cell_sets = list(T_cell      = c("CD3D", "CD3E", "CD3G", "CD2", "CD5", "LCK"),
                 CD8_T       = c("CD8A", "CD8B"),
                 Macrophage  = c("CD68", "CD163", "CD14", "C1QA", "C1QB"),
                 Endothelial = c("VWF", "PECAM1", "CDH5", "EMCN"),
                 Osteoclast  = c("ACP5", "CTSK", "MMP9"))
cell_sets = lapply(cell_sets, intersect, rownames(z_all))
lengths(cell_sets)
cs = sapply(cell_sets, function(g) colMeans(z_all[g, , drop = FALSE]))
sting = colMeans(z_all[sets21$T4_STING, ])

# 264. 各类细胞浸润评分：和 STING 的相关性，以及转移比较的 P 值
data.table(cell     = colnames(cs),
           rho_STING = round(apply(cs, 2, cor, y = sting, method = "spearman"), 2),
           p_meta    = signif(apply(cs, 2, p_meta), 3))

# 265. 逻辑回归：转移 ~ STING 评分 + T 细胞评分（看调整后 STING 的效应）
d = data.table(meta = as.integer(yes), STING = scale(sting)[, 1], Tcell = scale(cs[, "T_cell"])[, 1],
               Endo = scale(cs[, "Endothelial"])[, 1])
round(summary(glm(meta ~ STING,                 data = d, family = binomial))$coefficients, 3)
round(summary(glm(meta ~ STING + Tcell,         data = d, family = binomial))$coefficients, 3)
round(summary(glm(meta ~ STING + Tcell + Endo,  data = d, family = binomial))$coefficients, 3)
# 266. 同时放入巨噬细胞：看 STING 和巨噬细胞哪个更"领先"
round(summary(glm(meta ~ STING + Macro, data = cbind(d, Macro = scale(cs[, "Macrophage"])[, 1]),
                  family = binomial))$coefficients, 3)

# 267. 保存第十课结果
fwrite(data.table(sample = colnames(z_all), meta = yes, STING = sting, cs), "results/gse21257_immune_scores.csv")
# ===== 第十一课：FinnGen R13 重跑 MR =====

# 268. 下载 R13 的清单文件
dir.create("data_raw/finngen", showWarnings = FALSE)
mf_url  = "https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_manifest.tsv"
mf_file = "data_raw/finngen/finngen_R13_manifest.tsv"
system2("curl.exe", c("-L", "--ssl-no-revoke", "-o", mf_file, mf_url))
mf = fread(mf_file)
dim(mf); names(mf)

# 269. 找出和骨肉瘤、骨癌有关的表型，看病例数和对照数
hit = mf[grepl("OSTEOSARC|BONE", mf[[1]], ignore.case = TRUE) |
           grepl("osteosarcoma|bone and articular", mf[[2]], ignore.case = TRUE)]
hit[, 1:min(6, ncol(hit))]
# 270. 取出两个表型的下载网址
dl_tab = mf[phenocode %in% c("C3_OSTEOSARCOMA_WIDE", "C3_BONE_CARTILAGE"), .(phenocode, path_https)]
dl_tab

# 271. 逐个下载，断了自动续传，直到下完（可能要十几分钟到一小时）
for (i in seq_len(nrow(dl_tab))) {
  out = file.path("data_raw/finngen", basename(dl_tab$path_https[i]))
  repeat {
    st = system2("curl.exe", c("-L", "-C", "-", "--ssl-no-revoke", "--retry", "5", "-o", out, dl_tab$path_https[i]))
    if (st == 0) break
    Sys.sleep(5)
  }
}

# 272. 检查文件大小（单位 MB）
round(file.size(file.path("data_raw/finngen", basename(dl_tab$path_https))) / 1024^2, 1)
# 273. 加载包，重新准备 9 个锰 SNP 的暴露数据（和第 63 步一样）
library(data.table); library(TwoSampleMR)
iv  = fread("data_proc/mn_iv_main_clumped.csv")
exp = format_data(as.data.frame(iv), type = "exposure",
                  snp_col = "rsid", beta_col = "Effect", se_col = "StdErr",
                  effect_allele_col = "Allele1", other_allele_col = "Allele2",
                  eaf_col = "Freq1", pval_col = "P")
exp$exposure = "Manganese"; exp$id.exposure = "Mn"
nrow(exp)

# 274. 写一个函数：读 FinnGen 文件（只读需要的列），挑出这 9 个 SNP
read_fg = function(file, snps) {
  d = fread(file, select = c("rsids", "ref", "alt", "af_alt", "beta", "sebeta", "pval"))
  d[rsids %in% snps]
}
fg_os   = read_fg("data_raw/finngen/finngen_R13_C3_OSTEOSARCOMA_WIDE.gz", exp$SNP)
fg_bone = read_fg("data_raw/finngen/finngen_R13_C3_BONE_CARTILAGE.gz",    exp$SNP)
nrow(fg_os); nrow(fg_bone)
setdiff(exp$SNP, fg_os$rsids)

# 275. 整理成结局格式，统一方向（harmonisation）
fmt_out = function(d, name) {
  o = format_data(as.data.frame(d), type = "outcome", snp_col = "rsids",
                  beta_col = "beta", se_col = "sebeta", effect_allele_col = "alt",
                  other_allele_col = "ref", eaf_col = "af_alt", pval_col = "pval")
  o$outcome = name; o$id.outcome = name; o
}
dat_os   = harmonise_data(exp, fmt_out(fg_os,   "Osteosarcoma_R13"))
dat_bone = harmonise_data(exp, fmt_out(fg_bone, "BoneCancer_R13"))
table(dat_os$mr_keep); table(dat_bone$mr_keep)

# 276. 两个结局都跑三种 MR 方法，换算成 OR
mr_or = function(d) generate_odds_ratios(mr(d, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median")))[
  , c("outcome", "method", "nsnp", "or", "or_lci95", "or_uci95", "pval")]
rbind(mr_or(dat_os), mr_or(dat_bone))
# 277. 先确认 rs6541114 是否藏在"多个编号"的行里（rsids 一格写了好几个 rs 号）
fg_all = fread("data_raw/finngen/finngen_R13_C3_OSTEOSARCOMA_WIDE.gz", select = c("rsids", "ref", "alt", "af_alt", "beta", "sebeta", "pval"))
fg_all[grepl("rs6541114|rs10863510", rsids)]

# 278. 从完整的锰 GWAS 里找 rs10863510 的暴露数据（先查它的 GRCh37 位置）
library(ieugwasr)
pos_px = variants_rsid("rs10863510")
pos_px[, c("name", "chr", "pos")]
mn = readRDS("data_proc/mn_gwas_clean.rds")
px = mn[CHR == as.integer(pos_px$chr) & POS == as.integer(pos_px$pos)]
rm(mn); gc()
px[, .(MarkerName, Allele1, Allele2, Freq1, Effect, StdErr, P)]
# 274.（修正版）读取时去掉罕见变异（alt 频率 < 1% 或 > 99%），避免多等位基因位点混进来
read_fg = function(file, snps) {
  d = fread(file, select = c("rsids", "ref", "alt", "af_alt", "beta", "sebeta", "pval"))
  d[rsids %in% snps & af_alt > 0.01 & af_alt < 0.99]
}
rm(fg_all); gc()
# 279. 把 rs10863510 整理成暴露格式，替换掉 rs6541114
px[, rsid := "rs10863510"]
exp_px = format_data(as.data.frame(px), type = "exposure",
                     snp_col = "rsid", beta_col = "Effect", se_col = "StdErr",
                     effect_allele_col = "Allele1", other_allele_col = "Allele2",
                     eaf_col = "Freq1", pval_col = "P")
exp_px$exposure = "Manganese"; exp_px$id.exposure = "Mn"
exp9 = rbind(exp[exp$SNP != "rs6541114", ], exp_px[, names(exp)])
exp9$SNP

# 280. 重新从两个 FinnGen 文件里提取 9 个 SNP，统一方向
fg_os   = read_fg("data_raw/finngen/finngen_R13_C3_OSTEOSARCOMA_WIDE.gz", exp9$SNP)
fg_bone = read_fg("data_raw/finngen/finngen_R13_C3_BONE_CARTILAGE.gz",    exp9$SNP)
dat_os   = harmonise_data(exp9, fmt_out(fg_os,   "Osteosarcoma_R13"))
dat_bone = harmonise_data(exp9, fmt_out(fg_bone, "BoneCancer_R13"))
sum(dat_os$mr_keep); sum(dat_bone$mr_keep)

# 281. MR 主结果 + 异质性 + Egger 截距
rbind(mr_or(dat_os), mr_or(dat_bone))
mr_heterogeneity(dat_os)[, c("method", "Q", "Q_df", "Q_pval")]
mr_pleiotropy_test(dat_os)[, c("egger_intercept", "se", "pval")]

# 282. 保存
saveRDS(list(os = dat_os, bone = dat_bone), "data_proc/mr_dat_finngen_r13.rds")
# 283. 只用 2 个锰特异位点（SLC30A10 代替位点 rs10863510 + SLC39A8 rs13107325）
dat_os_mn = subset(dat_os, SNP %in% c("rs10863510", "rs13107325"))
generate_odds_ratios(mr(dat_os_mn, method_list = "mr_ivw"))[, c("method", "nsnp", "or", "or_lci95", "or_uci95", "pval")]

# 284. 每个 SNP 单独的估计（看哪些位点方向为正、哪些为负）
ss_os = mr_singlesnp(dat_os)
ss_os[, c("SNP", "b", "se", "p")]

# 285. 留一法：去掉任何一个 SNP 后结果是否改变
loo_os = mr_leaveoneout(dat_os)
loo_os[, c("SNP", "b", "se", "p")]

# 286. MR-PRESSO：检测离群 SNP
library(MRPRESSO)
set.seed(123)
presso_os = run_mr_presso(dat_os, NbDistribution = 3000)
presso_os[[1]]$`Main MR results`
presso_os[[1]]$`MR-PRESSO results`$`Global Test`$Pvalue

# 287. 画散点图和森林图，保存
p_sc = mr_scatter_plot(mr(dat_os, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median")), dat_os)[[1]]
p_fo = mr_forest_plot(ss_os)[[1]]
ggplot2::ggsave("results/mr_os_r13_scatter.pdf", p_sc, width = 6, height = 5)
ggplot2::ggsave("results/mr_os_r13_forest.pdf",  p_fo, width = 6, height = 5)
# ===== 第十二课：单细胞验证队列（GSE162454） =====

# 288. 列出 GEO 上 GSE162454 的补充文件
# [local-only, disabled for public release] setwd(".")
library(data.table)
v_dir = "data_raw/scRNA_val"
dir.create(v_dir, showWarnings = FALSE)
v_url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE162nnn/GSE162454/suppl/"
system2("curl.exe", c("-s", "--ssl-no-revoke", "-o", file.path(v_dir, "listing.html"), v_url))
lst = readLines(file.path(v_dir, "listing.html"))
gsub(".*href=\"([^\"]+)\".*", "\\1", grep("href=\"GSE", lst, value = TRUE))
p_sc
# 289. 诊断：看看下载到的网页有多大、写了什么
file.size(file.path(v_dir, "listing.html"))
length(lst)
head(lst, 20)

# 290. 直接查询 GSE162454_RAW.tar 是否存在（只看文件信息，不下载）
system2("curl.exe", c("-sI", "--ssl-no-revoke",
                      "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE162nnn/GSE162454/suppl/GSE162454_RAW.tar"))
# 291. 下载 GSE162454_RAW.tar（约 479 MB，断了自动续传）
v_tar = file.path(v_dir, "GSE162454_RAW.tar")
repeat {
  st = system2("curl.exe", c("-L", "-C", "-", "--ssl-no-revoke", "--retry", "5", "-o", v_tar,
                             "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE162nnn/GSE162454/suppl/GSE162454_RAW.tar"))
  if (st == 0) break
  Sys.sleep(5)
}
file.size(v_tar)

# 292. 解压，列出里面的文件
untar(v_tar, exdir = v_dir)
v_files = list.files(v_dir, recursive = TRUE)
v_files
# 297. 查询 GSE293065_RAW.tar 是否存在、有多大（只看信息，不下载）
sp_dir = "data_raw/spatial"
dir.create(sp_dir, showWarnings = FALSE)
system2("curl.exe", c("-sI", "--ssl-no-revoke",
                      "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE293nnn/GSE293065/suppl/GSE293065_RAW.tar"))
# 298. 下载 GSE293065_RAW.tar（约 396 MB，断了自动续传）
sp_tar = file.path(sp_dir, "GSE293065_RAW.tar")
repeat {
  st = system2("curl.exe", c("-L", "-C", "-", "--ssl-no-revoke", "--retry", "5", "-o", sp_tar,
                             "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE293nnn/GSE293065/suppl/GSE293065_RAW.tar"))
  if (st == 0) break
  Sys.sleep(5)
}
file.size(sp_tar)

# 299. 解压，列出里面的文件（只显示前 60 个）
untar(sp_tar, exdir = sp_dir)
sp_files = list.files(sp_dir, recursive = TRUE)
length(sp_files)
head(sp_files, 60)
# ===== 第十二课：单细胞验证（GSE162454）+ 空间数据准备 =====

# 300. 设置工作目录、加载包、定义文件夹
# [local-only, disabled for public release] setwd(".")
library(data.table); library(Seurat)
v_dir  = "data_raw/scRNA_val"
sp_dir = "data_raw/spatial"
file.size(file.path(sp_dir, "GSE293065_RAW.tar"))

# 299. 解压空间数据，看文件结构（只显示前 60 个）
untar(file.path(sp_dir, "GSE293065_RAW.tar"), exdir = sp_dir)
sp_files = list.files(sp_dir, recursive = TRUE)
length(sp_files)
head(sp_files, 60)

# 293. 找出单细胞验证队列 6 个样本的文件名前缀
v_pre = sub("_matrix.mtx.gz$", "", list.files(v_dir, pattern = "_matrix.mtx.gz$"))
v_pre

# 294. 逐个读入，建 Seurat 对象，合并
v_list = lapply(v_pre, function(p) {
  m = ReadMtx(mtx      = file.path(v_dir, paste0(p, "_matrix.mtx.gz")),
              cells    = file.path(v_dir, paste0(p, "_barcodes.tsv.gz")),
              features = file.path(v_dir, paste0(p, "_features.tsv.gz")))
  o = CreateSeuratObject(m, min.cells = 3, min.features = 200)
  o$patient = sub(".*_(OS_[0-9]+)$", "\\1", p)
  o
})
vs = merge(v_list[[1]], v_list[-1], add.cell.ids = sub(".*_(OS_[0-9]+)$", "\\1", v_pre))
rm(v_list); gc()
table(vs$patient)

# 295. 质控：和发现队列完全相同的标准
vs[["percent.mt"]] = PercentageFeatureSet(vs, pattern = "^MT-")
vq = as.data.table(vs@meta.data)
vq[, f_hi := 10^(median(log10(nFeature_RNA)) + 3 * mad(log10(nFeature_RNA))), by = patient]
vq[, keep := nFeature_RNA >= 300 & nFeature_RNA <= f_hi & percent.mt < 20]
vq[, .(f_hi = round(f_hi[1]), before = .N, after = sum(keep),
       kept = paste0(round(100 * mean(keep), 1), "%")), by = patient]
vs = vs[, vq$keep]
ncol(vs)
# 301. 去双细胞（按病人分别检测，要几分钟）
library(scDblFinder)
vs  = JoinLayers(vs)
sce = as.SingleCellExperiment(vs)
set.seed(123)
sce = scDblFinder(sce, samples = "patient")
vs$dbl = sce$scDblFinder.class
rm(sce); gc()
table(vs$patient, vs$dbl)
vs = subset(vs, dbl == "singlet")
ncol(vs)

# 302. 标准化 → 高变基因 → 缩放 → PCA
vs = NormalizeData(vs)
vs = FindVariableFeatures(vs, nfeatures = 2000)
vs = ScaleData(vs)
vs = RunPCA(vs, npcs = 30)

# 303. Harmony 按病人整合，然后聚类和 UMAP
library(harmony)
vs = RunHarmony(vs, group.by.vars = "patient")
vs = FindNeighbors(vs, reduction = "harmony", dims = 1:30)
vs = FindClusters(vs, resolution = 0.5)
vs = RunUMAP(vs, reduction = "harmony", dims = 1:30)

# 304. 用和发现队列完全相同的标志基因画气泡图
markers = list(
  Osteoblastic   = c("RUNX2", "ALPL", "IBSP", "COL1A1"),
  Chondroblastic = c("SOX9", "ACAN", "COL2A1"),
  Osteoclast     = c("ACP5", "CTSK", "MMP9"),
  Myeloid        = c("CD74", "LYZ", "CD68", "CD14"),
  T_NK           = c("CD3E", "IL7R", "NKG7", "GNLY"),
  B_Plasma       = c("MS4A1", "MZB1", "JCHAIN"),
  Fibroblast     = c("DCN", "LUM"),
  MSC            = c("CXCL12", "SFRP2", "MME"),
  Endothelial    = c("PECAM1", "VWF", "EGFL7"),
  Pericyte       = c("RGS5", "ACTA2"),
  Proliferating  = c("MKI67", "TOP2A"))
p_v = DotPlot(vs, features = markers, group.by = "seurat_clusters") + RotatedAxis()
ggplot2::ggsave("results/scRNA/val_dotplot_markers.pdf", p_v, width = 18, height = 7)
p_v

# 305. UMAP（按聚类、按病人），并保存对象
p1 = DimPlot(vs, group.by = "seurat_clusters", label = TRUE) + NoLegend()
p2 = DimPlot(vs, group.by = "patient")
p1 + p2
saveRDS(vs, "data_proc/sc_val_harmony.rds")
# 306. 每个群自动找最特异的基因（每群抽 500 个细胞，要几分钟）
set.seed(123)
vs_ds = subset(vs, downsample = 500)
Idents(vs_ds) = "seurat_clusters"
vmk   = FindAllMarkers(vs_ds, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
vtop5 = as.data.table(vmk)[order(cluster, -avg_log2FC)][, head(.SD, 5), by = cluster]
vtop5[, .(genes = paste(gene, collapse = ", ")), by = cluster]
rm(vs_ds); gc()
# 307. 把 26 个群对应到细胞类型
vanno = c("0" = "Myeloid", "3" = "Myeloid", "9" = "Myeloid", "14" = "Myeloid",
          "4" = "Osteoclast",
          "2" = "T_NK", "6" = "T_NK", "19" = "T_NK",
          "13" = "B_Plasma", "18" = "B_Plasma", "24" = "Mast",
          "10" = "OS_Chondroblastic",
          "11" = "Endothelial", "17" = "Pericyte", "23" = "Fibroblast_MSC")
os_cl = c("1", "5", "7", "8", "12", "15", "16", "20", "21", "22", "25")
vanno[os_cl] = "OS_Osteoblastic"
vs$celltype = unname(vanno[as.character(vs$seurat_clusters)])
table(vs$celltype, useNA = "ifany")

# 308. 带细胞类型标签的 UMAP，保存对象
p_vct = DimPlot(vs, group.by = "celltype", label = TRUE, repel = TRUE) + NoLegend()
ggplot2::ggsave("results/scRNA/val_umap_celltype.pdf", p_vct, width = 8, height = 7)
p_vct
saveRDS(vs, "data_proc/sc_val_annotated.rds")

# 309. 四层锰基因集（STING 基因的名字先确认是 TMEM173 还是 STING1）
c("TMEM173", "STING1", "ARG1") %in% rownames(vs)
sting_name = ifelse("STING1" %in% rownames(vs), "STING1", "TMEM173")
mn_sets = list(
  T1_Transport     = c("SLC39A8", "SLC39A14", "SLC30A10", "SLC11A2", "ATP2C1", "ATP13A1", "TMEM165", "TFRC"),
  T2_Enzyme        = c("SOD2", "ARG2", "GLUL", "PC", "PPM1A", "PPM1B"),
  T3_Glycosylation = c("B4GALT1", "GALNT1", "GALNT2", "XYLT1", "B3GAT3", "MGAT1"),
  T4_STING         = c("CGAS", sting_name, "TBK1", "IRF3", "IFNB1", "CXCL10", "ISG15"))
mn_sets = lapply(mn_sets, intersect, rownames(vs))
lengths(mn_sets)

# 310. 计算四层评分（和发现队列同样的方法）
vs = AddModuleScore(vs, features = mn_sets, name = "Tier")

# 311. 按病人取平均（伪批量），再做和第 185 步完全相同的 7 组配对比较
vmd = as.data.table(vs@meta.data)
vpb = vmd[, .(n = .N, T1_Transport = mean(Tier1), T2_Enzyme = mean(Tier2),
              T3_Glyco = mean(Tier3), T4_STING = mean(Tier4)), by = .(patient, celltype)][n >= 20]
vlong = melt(vpb, id.vars = c("patient", "celltype", "n"), variable.name = "tier", value.name = "score")
vwide = dcast(vlong, patient + tier ~ celltype, value.var = "score")
comps = data.table(
  tier = c("T1_Transport", "T1_Transport",      "T2_Enzyme", "T3_Glyco",       "T4_STING",    "T4_STING", "T4_STING"),
  cell = c("Osteoclast",   "OS_Chondroblastic", "Myeloid",   "Fibroblast_MSC", "Endothelial", "Pericyte", "T_NK"))
comps[, c("n_pt", "n_consistent", "median_diff", "p") := {
  d  = vwide[tier == .BY$tier]; x = d[[cell]]; y = d[["OS_Osteoblastic"]]
  ok = !is.na(x) & !is.na(y)
  list(sum(ok), sum(x[ok] > y[ok]), round(median(x[ok] - y[ok]), 3),
       if (sum(ok) >= 2) wilcox.test(x[ok], y[ok], paired = TRUE)$p.value else NA_real_)
}, by = .(tier, cell)]
comps[, p := signif(p, 3)]
comps
fwrite(comps, "results/scRNA/val_pseudobulk_paired_tests.csv")
# 312. 合并发现队列和验证队列：符号检验（方向一致的病人数 / 总病人数）
disc = fread("results/scRNA/pseudobulk_paired_tests.csv")[, .(tier, cell, n_disc = n_pt, k_disc = n_consistent)]
val  = fread("results/scRNA/val_pseudobulk_paired_tests.csv")[, .(tier, cell, n_val = n_pt, k_val = n_consistent)]
both = merge(disc, val, by = c("tier", "cell"))
both[, `:=`(n_all = n_disc + n_val, k_all = k_disc + k_val)]
both[, p_sign := signif(mapply(function(k, n) binom.test(k, n, 0.5)$p.value, k_all, n_all), 3)]
both[, p_sign_BH := signif(p.adjust(p_sign, "BH"), 3)]
both
fwrite(both, "results/scRNA/discovery_validation_combined.csv")
# ===== 第十三课：空间转录组（GSE293065） =====

# 313. 下载样本信息（series matrix），看 26 个样本各是什么（LR、MR 的含义）
sp_mx = file.path(sp_dir, "GSE293065_series_matrix.txt.gz")
system2("curl.exe", c("-L", "--ssl-no-revoke", "-o", sp_mx,
                      "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE293nnn/GSE293065/matrix/GSE293065_series_matrix.txt.gz"))
spl = readLines(gzfile(sp_mx))
grep("^!Sample_title|^!Sample_characteristics_ch1|^!Sample_source_name", spl, value = TRUE)

# 314. 安装读取需要的两个包（只需装一次）
install.packages(c("hdf5r", "R.utils"))

# 315. 把每个样本的文件整理成 Seurat 需要的文件夹结构（样本名/spatial/...），并解压
sp_pre = sub("_filtered_feature_bc_matrix.h5$", "", list.files(sp_dir, pattern = "_filtered_feature_bc_matrix.h5$"))
for (p in sp_pre) {
  d = file.path(sp_dir, "samples", p); dir.create(file.path(d, "spatial"), recursive = TRUE, showWarnings = FALSE)
  file.copy(file.path(sp_dir, paste0(p, "_filtered_feature_bc_matrix.h5")), file.path(d, "filtered_feature_bc_matrix.h5"))
  for (f in c("tissue_lowres_image.png", "tissue_hires_image.png", "scalefactors_json.json",
              "tissue_positions.csv", "tissue_positions_list.csv")) {
    src = file.path(sp_dir, paste0(p, "_", f, ".gz"))
    if (file.exists(src)) R.utils::gunzip(src, destname = file.path(d, "spatial", f), remove = FALSE, overwrite = TRUE)
  }
}
length(sp_pre)
list.files(file.path(sp_dir, "samples", sp_pre[1]), recursive = TRUE)

# 316. 试读第一个样本，画每个点的 UMI 数
library(Seurat)
st1 = Load10X_Spatial(file.path(sp_dir, "samples", sp_pre[1]))
dim(st1)
SpatialFeaturePlot(st1, features = "nCount_Spatial")
# 317. 整理样本信息：病灶类型、病人编号（MR10A/MR10B 合并为 MR10）
st_meta = data.table(title  = gsub('"', "", strsplit(grep("^!Sample_title", spl, value = TRUE), "\t")[[1]][-1]),
                     tissue = gsub('"|tissue: ', "", strsplit(grep("^!Sample_characteristics_ch1", spl, value = TRUE)[1], "\t")[[1]][-1]))
st_meta[, group   := fifelse(grepl("^PT", title), "Primary", fifelse(grepl("^LR", title), "Local_relapse", "Metastasis"))]
st_meta[, patient := sub("[AB]$", "", title)]
st_meta[, prefix  := sp_pre[match(title, sub("^GSM[0-9]+_", "", sp_pre))]]
st_meta[, .N, by = group]
uniqueN(st_meta$patient)

# 318. 定义细胞类型标志基因和锰四层基因（与单细胞相同；STING 用新名字）
cell_sets = list(OS_Osteoblastic = c("RUNX2", "ALPL", "IBSP", "COL1A1"),
                 Osteoclast      = c("ACP5", "CTSK", "MMP9", "DCSTAMP"),
                 Myeloid         = c("CD74", "LYZ", "CD68", "CD14"),
                 T_NK            = c("CD3E", "CD2", "NKG7", "GNLY"),
                 Fibroblast      = c("DCN", "LUM", "CXCL12", "SFRP2"),
                 Endothelial     = c("PECAM1", "VWF", "CDH5", "EGFL7"))
mn_sets = list(T1_Transport = c("SLC39A8", "SLC39A14", "SLC30A10", "SLC11A2", "ATP2C1", "ATP13A1", "TMEM165", "TFRC"),
               T2_Enzyme    = c("SOD2", "ARG2", "GLUL", "PC", "PPM1A", "PPM1B"),
               T3_Glyco     = c("B4GALT1", "GALNT1", "GALNT2", "XYLT1", "B3GAT3", "MGAT1"),
               T4_STING     = c("CGAS", "STING1", "TMEM173", "TBK1", "IRF3", "IFNB1", "CXCL10", "ISG15"))
pairs = data.table(cell = c("Osteoclast",   "Osteoclast", "Myeloid",   "Fibroblast", "T_NK",     "Endothelial"),
                   mn   = c("T1_Transport", "aVb3",       "T2_Enzyme", "T3_Glyco",   "T4_STING", "T4_STING"))

# 319. 每个样本：读入 → 标准化 → 打分 → 算预设的 6 对相关（要 10–20 分钟）
score_one = function(p) {
  st   = Load10X_Spatial(file.path(sp_dir, "samples", p))
  st   = NormalizeData(st)
  sets = lapply(c(cell_sets, mn_sets), intersect, rownames(st))
  st   = AddModuleScore(st, features = sets, name = "S")
  sc   = as.data.table(st@meta.data[, paste0("S", seq_along(sets))]); setnames(sc, names(sets))
  ex   = GetAssayData(st, layer = "data")
  sc[, aVb3 := colMeans(as.matrix(ex[intersect(c("ITGAV", "ITGB3"), rownames(ex)), , drop = FALSE]))]
  out  = copy(pairs)
  out[, rho := mapply(function(a, b) cor(sc[[a]], sc[[b]], method = "spearman"), cell, mn)]
  out[, `:=`(prefix = p, n_spots = ncol(st))]
  out
}
sp_cor = rbindlist(lapply(st_meta$prefix, score_one))
sp_cor = merge(sp_cor, st_meta[, .(prefix, title, patient, group)], by = "prefix")
fwrite(sp_cor, "results/spatial_within_sample_cor.csv")

# 320. 按病人汇总（同一病人的 A/B 取平均），看 24 个病人里有几个是正相关
sp_pt = sp_cor[, .(rho = mean(rho)), by = .(patient, cell, mn)]
sp_pt[, .(n_pt = .N, n_pos = sum(rho > 0), median_rho = round(median(rho), 2),
          p_sign = signif(binom.test(sum(rho > 0), .N, 0.5)$p.value, 3)), by = .(cell, mn)]
sapply(c(cell_sets, mn_sets), function(g) length(intersect(g, rownames(st1))))
# 322. 重跑：每个样本算所有"细胞类型 × 锰层级"的相关（用于特异性对照）
all_pairs = CJ(cell = names(cell_sets), mn = c(names(mn_sets), "aVb3"))
pairs = all_pairs
sp_all = rbindlist(lapply(st_meta$prefix, score_one))
sp_all = merge(sp_all, st_meta[, .(prefix, patient, group)], by = "prefix")
fwrite(sp_all, "results/spatial_all_pairs_cor.csv")

# 323. 按病人平均后，做成"细胞 × 层级"的中位 ρ 矩阵，画热图
sp_allpt = sp_all[, .(rho = mean(rho)), by = .(patient, cell, mn)]
rho_mat  = dcast(sp_allpt[, .(rho = median(rho)), by = .(cell, mn)], cell ~ mn, value.var = "rho")
rm_mat   = as.matrix(rho_mat[, -1]); rownames(rm_mat) = rho_mat$cell
round(rm_mat, 2)
pheatmap::pheatmap(rm_mat, cluster_rows = FALSE, cluster_cols = FALSE,
                   display_numbers = TRUE, number_format = "%.2f",
                   filename = "results/spatial_rho_heatmap.pdf", width = 6, height = 4.5)
# 324. 出版版热图：0 为白色，正负对称
lim = max(abs(rm_mat))
brk = seq(-lim, lim, length.out = 101)
pal = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
rm_ord = rm_mat[c("OS_Osteoblastic", "Osteoclast", "Myeloid", "Fibroblast", "Endothelial", "T_NK"),
                c("T1_Transport", "T2_Enzyme", "T3_Glyco", "T4_STING", "aVb3")]
pheatmap::pheatmap(rm_ord, cluster_rows = FALSE, cluster_cols = FALSE, color = pal, breaks = brk,
                   display_numbers = TRUE, number_format = "%.2f", border_color = "grey80",
                   main = "Median within-sample Spearman rho (24 patients)",
                   filename = "results/spatial_rho_heatmap_pub.pdf", width = 6, height = 4.5)
# 324.（修正）两端颜色稍浅，让数字看得清；加宽避免标题被截断
pal = colorRampPalette(c("#4393C3", "white", "#D6604D"))(100)
pheatmap::pheatmap(rm_ord, cluster_rows = FALSE, cluster_cols = FALSE, color = pal, breaks = brk,
                   display_numbers = TRUE, number_format = "%.2f", fontsize_number = 9,
                   border_color = "grey80", main = "Median within-sample Spearman rho (24 patients)",
                   filename = "results/spatial_rho_heatmap_pub.pdf", width = 7, height = 4.5)
# 325. 每个样本：选出肿瘤富集点（成骨型评分前 25%），计算四层基因在这些点里的平均表达
tumor_one = function(p) {
  st   = NormalizeData(Load10X_Spatial(file.path(sp_dir, "samples", p)))
  ex   = GetAssayData(st, layer = "data")
  os_g = intersect(cell_sets$OS_Osteoblastic, rownames(ex))
  os_s = colMeans(as.matrix(ex[os_g, , drop = FALSE]))
  tum  = os_s >= quantile(os_s, 0.75)
  res  = sapply(mn_sets, function(g) { g = intersect(g, rownames(ex)); mean(as.matrix(ex[g, tum, drop = FALSE])) })
  data.table(prefix = p, n_tumor_spots = sum(tum), t(res))
}
sp_tum = rbindlist(lapply(st_meta$prefix, tumor_one))
sp_tum = merge(sp_tum, st_meta[, .(prefix, patient, group)], by = "prefix")

# 326. 按病人合并（A/B 取平均），分成原发 vs 复发/转移
tum_pt = sp_tum[, lapply(.SD, mean), by = .(patient, group), .SDcols = names(mn_sets)]
tum_pt[, grp2 := fifelse(group == "Primary", "Primary", "Relapse_Met")]
tum_pt[, .N, by = grp2]
fwrite(tum_pt, "results/spatial_tumor_spots_by_patient.csv")

# 327. 比较：T1 是预设的主要假设，其他三层为探索性
tl = melt(tum_pt, id.vars = c("patient", "grp2"), measure.vars = names(mn_sets),
          variable.name = "tier", value.name = "expr")
tl[, .(Primary     = round(median(expr[grp2 == "Primary"]), 3),
       Relapse_Met = round(median(expr[grp2 == "Relapse_Met"]), 3),
       p = signif(wilcox.test(expr[grp2 == "Relapse_Met"], expr[grp2 == "Primary"])$p.value, 3)), by = tier]

# 328. 画图：每个点一个病人
p_tum = ggplot2::ggplot(tl, ggplot2::aes(grp2, expr, color = grp2)) +
  ggplot2::geom_boxplot(outlier.shape = NA) + ggplot2::geom_jitter(width = 0.15, size = 2) +
  ggplot2::facet_wrap(~ tier, scales = "free_y", ncol = 4) + ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "none") + ggplot2::labs(x = NULL, y = "Mean expression in tumour-rich spots")
ggplot2::ggsave("results/spatial_tumor_primary_vs_met.pdf", p_tum, width = 11, height = 4)
p_tum
# 329. 选代表性样本：髓系–锰依赖酶的相关系数最接近中位数、点数 > 1000
rep_cand = sp_all[cell == "Myeloid" & mn == "T2_Enzyme" & n_spots > 1000]
rep_cand[, dist := abs(rho - median(rho))]
rep_p = rep_cand[order(dist)][1, prefix]
rep_cand[order(dist)][1:3, .(prefix, rho, n_spots)]

# 330. 读入这个样本，计算要展示的几个评分
st = NormalizeData(Load10X_Spatial(file.path(sp_dir, "samples", rep_p)))
show_sets = list(Myeloid = cell_sets$Myeloid, T2_Enzyme = mn_sets$T2_Enzyme,
                 Osteoclast = cell_sets$Osteoclast, Tumour = cell_sets$OS_Osteoblastic)
show_sets = lapply(show_sets, intersect, rownames(st))
st = AddModuleScore(st, features = show_sets, name = "Show")
names(st@meta.data)[grep("^Show", names(st@meta.data))] = names(show_sets)
ex = GetAssayData(st, layer = "data")
st$aVb3 = colMeans(as.matrix(ex[intersect(c("ITGAV", "ITGB3"), rownames(ex)), , drop = FALSE]))

# 331. 画 6 张并排的空间图：H&E、肿瘤、髓系、锰依赖酶、破骨、αvβ3
p_he = SpatialFeaturePlot(st, features = "nCount_Spatial", alpha = 0) + ggplot2::ggtitle("H&E")
p_sp = SpatialFeaturePlot(st, features = c("Tumour", "Myeloid", "T2_Enzyme", "Osteoclast", "aVb3"),
                          ncol = 5, pt.size.factor = 1.6)
p_he
p_sp
ggplot2::ggsave("results/spatial_example_maps.pdf", p_sp, width = 20, height = 5)
p_tum
# 332. 主图：24 个病人各自的相关系数（6 个预设假设）
plot_pt = sp_pt[, .(patient, pair = paste(cell, "~", mn), rho)]
p_rho = ggplot2::ggplot(plot_pt, ggplot2::aes(x = rho, y = pair)) +
  ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  ggplot2::geom_jitter(height = 0.15, size = 1.8, alpha = 0.7, colour = "#2166AC") +
  ggplot2::stat_summary(fun = median, geom = "point", shape = 124, size = 8, colour = "#B2182B") +
  ggplot2::theme_bw() + ggplot2::labs(x = "Within-sample Spearman rho (per patient)", y = NULL)
p_rho
ggplot2::ggsave("results/spatial_rho_per_patient.pdf", p_rho, width = 7, height = 3.5)

# 333. 邻域平滑：每个点取自己和最近 6 个点的平均（只用于作图）
coords = GetTissueCoordinates(st)
coords = coords[colnames(st), ]
D  = as.matrix(dist(coords[, c("x", "y")]))
nn = t(apply(D, 1, function(r) order(r)[1:7]))
rm(D); gc()
smooth = function(v) rowMeans(matrix(v[nn], nrow = nrow(nn)))
for (f in c("Tumour", "Myeloid", "T2_Enzyme", "Osteoclast", "aVb3")) st[[paste0(f, "_s")]] = smooth(st[[f]][, 1])

# 334. 画平滑后的空间图
p_sp2 = SpatialFeaturePlot(st, features = paste0(c("Tumour", "Myeloid", "T2_Enzyme", "Osteoclast", "aVb3"), "_s"),
                           ncol = 5, pt.size.factor = 1.6)
p_sp2
ggplot2::ggsave("results/spatial_example_maps_smoothed.pdf", p_sp2, width = 20, height = 5)
file.exists(c("results/spatial_rho_per_patient.pdf", "results/spatial_example_maps_smoothed.pdf"))
graphics.off()
print(p_rho)
print(p_sp2)
# ===== 第十四课：药物靶点 =====

# 335. 预设靶点清单（全部来自本研究的发现）
targets = data.table(
  gene = c("SOD2", "GLUL", "ARG2",
           "CGAS", "STING1", "TBK1", "IRF3",
           "ITGAV", "ITGB3", "ITGB1", "ITGA10",
           "MIF", "CD74", "PPIA", "BSG", "SPP1", "CD44", "SEMA4D", "PLXNB1",
           "SLC39A8", "SLC39A14", "SLC11A2", "TFRC",
           "NOD2", "CSF1R"),
  module = c(rep("Mn_enzyme_myeloid", 3), rep("STING_pathway", 4), rep("Integrin", 4),
             rep("Cell_communication", 8), rep("Mn_transporter", 4), rep("Macrophage_activation", 2)))
targets[, .N, by = module]

# 336. 写好 DGIdb 的查询内容（GraphQL 格式），存成文件
dir.create("data_raw/drug", showWarnings = FALSE)
q = paste0('{"query":"{ genes(names: [', paste0('\\"', targets$gene, '\\"', collapse = ","),
           ']) { nodes { name interactions { drug { name approved } interactionScore ',
           'interactionTypes { type } sources { sourceDbName } } } } }"}')
writeLines(q, "data_raw/drug/dgidb_query.json")

# 337. 发送查询，保存结果
system2("curl.exe", c("-s", "--ssl-no-revoke", "-X", "POST",
                      "-H", shQuote("Content-Type: application/json", type = "cmd"),
                      "--data", "@data_raw/drug/dgidb_query.json",
                      "-o", "data_raw/drug/dgidb_result.json", "https://dgidb.org/api/graphql"))
file.size("data_raw/drug/dgidb_result.json")

# 338. 整理成表格：每行一个"基因–药物"组合
res   = jsonlite::fromJSON("data_raw/drug/dgidb_result.json", simplifyVector = FALSE)
nodes = res$data$genes$nodes
dgi = rbindlist(lapply(nodes, function(g) {
  if (length(g$interactions) == 0) return(NULL)
  rbindlist(lapply(g$interactions, function(it) data.table(
    gene     = g$name,
    drug     = it$drug$name,
    approved = isTRUE(it$drug$approved),
    score    = if (is.null(it$interactionScore)) NA_real_ else it$interactionScore,
    type     = paste(unique(sapply(it$interactionTypes, `[[`, "type")), collapse = ";"),
    n_source = length(it$sources))))
}))
dgi = merge(dgi, targets, by = "gene")
dgi[, .(n_drugs = uniqueN(drug), n_approved = uniqueN(drug[approved])), by = .(module, gene)][order(module, -n_drugs)]
fwrite(dgi, "results/drug/dgidb_interactions.csv")
dir.create("results/drug", showWarnings = FALSE)
fwrite(dgi, "results/drug/dgidb_interactions.csv")
# 339. 只看标明了作用类型的组合，每个靶点按相互作用评分取前 3 个药物
dgi_typed = dgi[type != ""]
dgi_typed[order(gene, -score)][, head(.SD, 3), by = gene][, .(module, gene, drug, type, approved, score = round(score, 2))]

# 340. 专门查和我们发现直接相关的几类药物
key = "MIFAMURTIDE|CILENGITIDE|MANGANESE|ADU-S100|MIW815|DIABZI|PEXIDARTINIB|DENOSUMAB|ZOLEDRON"
dgi[grepl(key, toupper(drug)), .(gene, drug, type, approved, n_source)]
# 341. 整理一张精简的"关键靶点–药物"补充表
key_tab = dgi[grepl(key, toupper(drug)) | (type != "" & gene %in% c("ARG2", "CD74", "MIF", "BSG", "SEMA4D", "TBK1", "CGAS")),
              .(drugs = paste(unique(drug), collapse = "; "),
                types = paste(unique(type[type != ""]), collapse = "; "),
                any_approved = any(approved)), by = .(module, gene)]
key_tab
fwrite(key_tab, "results/drug/key_target_drug_table.csv")

# 342. 发现队列：25 个靶点在各类细胞中的表达（STING1 在这批数据里叫 TMEM173）
sc_tme  = readRDS("data_proc/sc_tme_scored.rds")
tg_disc = sub("^STING1$", "TMEM173", targets$gene)
tg_disc = intersect(tg_disc, rownames(sc_tme))
p_t1 = DotPlot(sc_tme, features = tg_disc, group.by = "celltype") + RotatedAxis() + ggplot2::ggtitle("Discovery (GSE152048)")
ggplot2::ggsave("results/drug/target_dotplot_discovery.pdf", p_t1, width = 12, height = 5)
rm(sc_tme); gc()

# 343. 验证队列：同样的靶点
vs = readRDS("data_proc/sc_val_annotated.rds")
tg_val = intersect(sub("^STING1$", ifelse("STING1" %in% rownames(vs), "STING1", "TMEM173"), targets$gene), rownames(vs))
p_t2 = DotPlot(vs, features = tg_val, group.by = "celltype") + RotatedAxis() + ggplot2::ggtitle("Validation (GSE162454)")
ggplot2::ggsave("results/drug/target_dotplot_validation.pdf", p_t2, width = 12, height = 5)
print(p_t1)
print(p_t2)
# 344. 在 ClinicalTrials.gov 查骨肉瘤相关的试验：按药物/靶点分别检索
ct_q = c(mifamurtide = "mifamurtide OR MTP-PE", pexidartinib = "pexidartinib OR CSF1R",
         cilengitide = "cilengitide OR integrin", sting = "STING agonist",
         manganese = "manganese", zoledronic = "zoledronic acid", denosumab = "denosumab")
ct = rbindlist(lapply(names(ct_q), function(k) {
  url = paste0("https://clinicaltrials.gov/api/v2/studies?format=json&pageSize=100",
               "&query.cond=osteosarcoma&query.intr=", URLencode(ct_q[[k]], reserved = TRUE),
               "&fields=NCTId,BriefTitle,Phase,OverallStatus,StartDate")
  f = file.path("data_raw/drug", paste0("ct_", k, ".json"))
  system2("curl.exe", c("-s", "--ssl-no-revoke", "-o", f, shQuote(url, type = "cmd")))
  s = jsonlite::fromJSON(f)$studies
  if (length(s) == 0) return(data.table(query = k, n = 0))
  data.table(query  = k,
             nct    = s$protocolSection$identificationModule$nctId,
             title  = s$protocolSection$identificationModule$briefTitle,
             phase  = sapply(s$protocolSection$designModule$phases, paste, collapse = "/"),
             status = s$protocolSection$statusModule$overallStatus)
}), fill = TRUE)

# 345. 每类药物有几个骨肉瘤试验
ct[!is.na(nct), .N, by = query]
fwrite(ct, "results/drug/clinicaltrials_osteosarcoma.csv")
# 346. 查看每个试验的标题、分期、状态
ct[!is.na(nct), .(query, nct, phase, status, title)]
# 347. 人工标注相关性，保存最终版本
rel = c(NCT03737435 = "Yes", NCT00631631 = "Yes", NCT04383288 = "Check", NCT07787429 = "Yes",
        NCT03643133 = "Yes", NCT02441309 = "Yes", NCT01194284 = "Check", NCT04571229 = "Yes",
        NCT01459484 = "Check", NCT05093322 = "Yes_CSF1R_surufatinib",
        NCT02517918 = "Check", NCT05400603 = "No_gdT_expansion", NCT00691236 = "Yes", NCT03932071 = "Yes",
        NCT00742924 = "Yes", NCT00470223 = "Yes", NCT02508038 = "No_transplant",
        NCT00987636 = "No_Ewing", NCT06625190 = "No_transplant", NCT02470091 = "Yes")
ct[, relevant := rel[nct]]
ct[!is.na(nct), .N, by = .(query, relevant)][order(query)]
fwrite(ct, "results/drug/clinicaltrials_osteosarcoma.csv")
# 348. 检查已有图片文件是否都存在
figs = c("results/mr_os_r13_scatter.pdf", "results/mr_os_r13_forest.pdf",
         "results/mr_bone_scatter.pdf", "results/mr_bone_forest.pdf",
         "results/scRNA/umap_celltype.pdf", "results/scRNA/dotplot_mn_genes.pdf",
         "results/scRNA/heatmap_mn_tiers.pdf", "results/scRNA/pseudobulk_tiers_boxplot.pdf",
         "results/scRNA/val_umap_celltype.pdf", "results/scRNA/qc_violin_before.pdf",
         "results/scRNA/dotplot_markers.pdf", "results/scRNA/val_dotplot_markers.pdf",
         "results/scRNA/os_lesion_compare.pdf",
         "results/scRNA/cellchat_circle.pdf", "results/scRNA/cellchat_integrin_heatmap.pdf",
         "results/scRNA/cellchat_bubble_mn_to_os.pdf", "results/scRNA/dotplot_receptors.pdf",
         "results/spatial_rho_heatmap_pub.pdf", "results/spatial_rho_per_patient.pdf",
         "results/spatial_example_maps_smoothed.pdf", "results/spatial_tumor_primary_vs_met.pdf",
         "results/drug/target_dotplot_discovery.pdf", "results/drug/target_dotplot_validation.pdf")
data.table(file = figs, exists = file.exists(figs), KB = round(file.size(figs) / 1024))
# 350.（修正）Fig 1B：只保留 R13，按结局分组，数字放在图外的一列
fig1b[, method := sub("Inverse variance weighted", "IVW", method)]
fig1b[, method := sub("MR Egger", "MR-Egger", method)]
fig1b[, group := fifelse(grepl("^Bone", outcome), "Bone cancer", "Osteosarcoma")]
main_b = fig1b[grepl("R13|Mn-specific", outcome) & !(grepl("^Bone", outcome) & method != "IVW")]
main_b[, method := fifelse(grepl("Mn-specific", outcome), "IVW, Mn-specific SNPs (SLC30A10, SLC39A8)", method)]
main_b[, group := fifelse(group == "Osteosarcoma", "Osteosarcoma\n(FinnGen R13, 134 cases)", "Bone cancer\n(FinnGen R13, 338 cases)")]
main_b[, group := factor(group, levels = unique(group))]
main_b[, method := factor(method, levels = rev(unique(method)))]

p1b = ggplot(main_b, aes(x = or, y = method)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_pointrange(aes(xmin = lo, xmax = hi), size = 0.4, colour = "#2166AC") +
  geom_text(aes(x = 6, label = label), hjust = 0, size = 3.2) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4), limits = c(0.25, 4)) +
  coord_cartesian(clip = "off") +
  theme_bw() +
  theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0),
        plot.margin = margin(5, 110, 5, 5)) +
  labs(x = "OR per SD increase in blood Mn (log scale)", y = NULL)
p1b
ggsave("results/figures/Fig1B_MR_forest.pdf", p1b, width = 8, height = 3.5)

# 350b. 补充图：R5 骨癌（旧数据）
sup_b = fig1b[grepl("R5", outcome)]
pS = ggplot(sup_b, aes(x = or, y = method)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_pointrange(aes(xmin = lo, xmax = hi), size = 0.4) +
  geom_text(aes(x = 8, label = label), hjust = 0, size = 3.2) +
  scale_x_log10(breaks = c(0.5, 1, 2, 4)) + coord_cartesian(xlim = c(0.4, 6), clip = "off") +
  theme_bw() + theme(plot.margin = margin(5, 110, 5, 5)) +
  labs(x = "OR (log scale)", y = NULL, title = "Bone cancer, FinnGen R5 (119 cases)")
ggsave("results/figures/FigS_MR_R5.pdf", pS, width = 7, height = 2.5)

1
graphics.off()
plot(1:10)
# 350.（完整版）整理数据 → 画 Fig 1B → 保存 → 显示
fig1b[, method := sub("Inverse variance weighted", "IVW", method)]
fig1b[, method := sub("MR Egger", "MR-Egger", method)]

main_b = fig1b[(grepl("^Osteosarcoma", outcome)) |
                 (grepl("^Bone cancer, FinnGen R13", outcome) & method == "IVW")]
main_b[, method := as.character(method)]
main_b[grepl("Mn-specific", outcome), method := "IVW, Mn-specific SNPs"]
main_b[, group := fifelse(grepl("^Osteosarcoma", outcome),
                          "Osteosarcoma\n(R13, 134 cases)", "Bone cancer\n(R13, 338 cases)")]
main_b[, group  := factor(group, levels = c("Osteosarcoma\n(R13, 134 cases)", "Bone cancer\n(R13, 338 cases)"))]
main_b[, method := factor(method, levels = rev(c("IVW", "MR-Egger", "Weighted median", "IVW, Mn-specific SNPs")))]
main_b[, .(group, method, label)]

p1b = ggplot(main_b, aes(x = or, y = method)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_pointrange(aes(xmin = lo, xmax = hi), size = 0.4, colour = "#2166AC") +
  geom_text(aes(x = 6, label = label), hjust = 0, size = 3.2) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4)) +
  coord_cartesian(xlim = c(0.25, 4), clip = "off") +
  theme_bw() +
  theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0),
        plot.margin = margin(5, 110, 5, 5)) +
  labs(x = "OR per SD increase in blood Mn (log scale)", y = NULL)
ggsave("results/figures/Fig1B_MR_forest.pdf", p1b, width = 8, height = 3.5)
print(p1b)
# 350c.（可选美化）分组标签去掉灰底，数字列加表头
p1b = p1b +
  theme(strip.background = element_blank(), strip.text.y.left = element_text(angle = 0, face = "bold")) +
  annotate("text", x = 6, y = Inf, label = "OR (95% CI)", hjust = 0, vjust = -0.6, size = 3.2, fontface = "bold")
ggsave("results/figures/Fig1B_MR_forest.pdf", p1b, width = 8, height = 3.5)
print(p1b)
# 350c.（修正）表头只在骨肉瘤组上方显示一次
hdr = data.frame(group = factor(levels(main_b$group)[1], levels = levels(main_b$group)),
                 x = 6, y = Inf, label = "OR (95% CI)")
p1b = ggplot(main_b, aes(x = or, y = method)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_pointrange(aes(xmin = lo, xmax = hi), size = 0.4, colour = "#2166AC") +
  geom_text(aes(x = 6, label = label), hjust = 0, size = 3.2) +
  geom_text(data = hdr, aes(x = x, y = y, label = label), inherit.aes = FALSE,
            hjust = 0, vjust = -0.8, size = 3.2, fontface = "bold") +
  facet_grid(group ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4)) +
  coord_cartesian(xlim = c(0.25, 4), clip = "off") +
  theme_bw() +
  theme(strip.placement = "outside", strip.background = element_blank(),
        strip.text.y.left = element_text(angle = 0, face = "bold"),
        plot.margin = margin(20, 110, 5, 5)) +
  labs(x = "OR per SD increase in blood Mn (log scale)", y = NULL)
ggsave("results/figures/Fig1B_MR_forest.pdf", p1b, width = 8, height = 3.5)
print(p1b)
# 351. Fig 3C：两个单细胞队列中，每个比较"方向一致的病人比例"
cv = fread("results/scRNA/discovery_validation_combined.csv")
cv[, comp := paste0(cell, " · ", tier)]
cv_long = rbind(cv[, .(comp, cohort = "Discovery (n=11)", k = k_disc, n = n_disc)],
                cv[, .(comp, cohort = "Validation (n=6)", k = k_val,  n = n_val)])
cv_long[, prop := k / n]
cv[, star := sprintf("%d/%d, BH P = %.1e", k_all, n_all, p_sign_BH)]
p3c = ggplot(cv_long, aes(x = prop, y = comp, fill = cohort)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  geom_text(aes(label = paste0(k, "/", n)), position = position_dodge(0.8), hjust = -0.1, size = 3) +
  geom_text(data = cv, aes(x = 1.35, y = comp, label = star), inherit.aes = FALSE, hjust = 1, size = 3) +
  geom_vline(xintercept = 0.5, linetype = 2, colour = "grey50") +
  scale_x_continuous(limits = c(0, 1.35), breaks = c(0, 0.5, 1), labels = scales::percent) +
  scale_fill_manual(values = c("#4393C3", "#D6604D")) +
  theme_bw() + labs(x = "Patients with concordant direction (vs osteoblastic OS cells)", y = NULL, fill = NULL)
ggsave("results/figures/Fig3C_concordance.pdf", p3c, width = 10, height = 4)
print(p3c)
# 351.（修正）Fig 3C：合并统计单独放一列，标签可读，按证据强度排序
lab_map = c("Myeloid · T2_Enzyme"           = "Myeloid — Mn enzymes",
            "Endothelial · T4_STING"        = "Endothelial — STING",
            "Osteoclast · T1_Transport"     = "Osteoclast — Mn transport",
            "Fibroblast_MSC · T3_Glyco"     = "Fibroblast/MSC — Mn glycosylation",
            "Pericyte · T4_STING"           = "Pericyte — STING",
            "T_NK · T4_STING"               = "T/NK — STING",
            "OS_Chondroblastic · T1_Transport" = "Chondroblastic OS — Mn transport")
cv[, comp2 := factor(lab_map[comp], levels = rev(lab_map))]
cv_long[, comp2 := factor(lab_map[comp], levels = rev(lab_map))]
cv_long[, cohort := factor(cohort, levels = c("Validation (n=6)", "Discovery (n=11)"))]
cv[, star := sprintf("%d/%d   %s", k_all, n_all,
                     ifelse(p_sign_BH < 0.001, formatC(p_sign_BH, format = "e", digits = 1),
                            formatC(p_sign_BH, format = "f", digits = 3)))]

p3c = ggplot(cv_long, aes(x = prop, y = comp2, fill = cohort)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  geom_text(aes(label = paste0(k, "/", n)), position = position_dodge(0.8), hjust = 1.15,
            size = 2.8, colour = "white") +
  geom_text(data = cv, aes(x = 1.06, y = comp2, label = star), inherit.aes = FALSE, hjust = 0, size = 3) +
  annotate("text", x = 1.06, y = Inf, label = "Pooled         BH P", hjust = 0, vjust = -0.8,
           size = 3, fontface = "bold") +
  geom_vline(xintercept = 0.5, linetype = 2, colour = "grey50") +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = scales::percent) +
  coord_cartesian(xlim = c(0, 1), clip = "off") +
  scale_fill_manual(values = c("Discovery (n=11)" = "#4393C3", "Validation (n=6)" = "#D6604D"),
                    breaks = c("Discovery (n=11)", "Validation (n=6)")) +
  theme_bw() +
  theme(legend.position = "bottom", plot.margin = margin(20, 150, 5, 5)) +
  labs(x = "Patients with concordant direction (vs osteoblastic OS cells)", y = NULL, fill = NULL)
ggsave("results/figures/Fig3C_concordance.pdf", p3c, width = 9, height = 4.5)
print(p3c)
p3c$layers[[4]]$aes_params$label = "Pooled         BH P"
ggsave("results/figures/Fig3C_concordance.pdf", p3c, width = 9, height = 4.5)
# 352. 两个队列的细胞浸润评分（同一套标志基因、同样的 z 分数平均法）
cell_sets_b = list(T_cell      = c("CD3D", "CD3E", "CD3G", "CD2", "CD5", "LCK"),
                   CD8_T       = c("CD8A", "CD8B"),
                   Macrophage  = c("CD68", "CD163", "CD14", "C1QA", "C1QB"),
                   Endothelial = c("VWF", "PECAM1", "CDH5", "EMCN"),
                   Osteoclast  = c("ACP5", "CTSK", "MMP9"))
mat = readRDS("data_proc/target_expr.rds")
bk  = fread("data_proc/target_clin_scores.csv")
zt  = t(scale(t(mat[intersect(unlist(cell_sets_b), rownames(mat)), bk$sample])))
cs_t = sapply(cell_sets_b, function(g) colMeans(zt[intersect(g, rownames(zt)), , drop = FALSE]))
g21i = fread("results/gse21257_immune_scores.csv")

# 353. Fig 6A：STING 评分与各类细胞浸润的 Spearman 相关（两队列并排）
rho6 = rbind(data.table(cohort = "TARGET (RNA-seq, n=88)", cell = names(cell_sets_b),
                        rho = apply(cs_t, 2, cor, y = bk$T4_STING, method = "spearman")),
             data.table(cohort = "GSE21257 (array, n=53)", cell = names(cell_sets_b),
                        rho = sapply(names(cell_sets_b), function(k) cor(g21i[[k]], g21i$STING, method = "spearman"))))
rho6[, cell := factor(cell, levels = rev(names(cell_sets_b)))]
p6a = ggplot(rho6, aes(x = cohort, y = cell, fill = rho)) +
  geom_tile(colour = "white") + geom_text(aes(label = sprintf("%.2f", rho)), size = 3.5) +
  scale_fill_gradient2(low = "#4393C3", mid = "white", high = "#D6604D", limits = c(-1, 1)) +
  theme_minimal() + labs(x = NULL, y = NULL, fill = "Spearman\nrho", title = "STING score vs cell infiltration")
print(p6a)
ggsave("results/figures/Fig6A_STING_cell_cor.pdf", p6a, width = 5, height = 4)

# 354. Fig 6B：GSE21257 逻辑回归，STING 在调整 T 细胞 / 巨噬细胞前后的 OR
g21i[, `:=`(meta01 = as.integer(meta), zS = scale(STING)[, 1], zT = scale(T_cell)[, 1], zM = scale(Macrophage)[, 1])]
get_or = function(f, model) {
  s = summary(glm(f, data = g21i, family = binomial))$coefficients[-1, , drop = FALSE]
  data.table(model, term = rownames(s), or = exp(s[, 1]),
             lo = exp(s[, 1] - 1.96 * s[, 2]), hi = exp(s[, 1] + 1.96 * s[, 2]), p = s[, 4])
}
or6 = rbind(get_or(meta01 ~ zS,      "STING alone"),
            get_or(meta01 ~ zS + zT, "STING + T cells"),
            get_or(meta01 ~ zS + zM, "STING + macrophages"))
or6[, term := c(zS = "STING score", zT = "T-cell score", zM = "Macrophage score")[term]]
or6[, model := factor(model, levels = rev(c("STING alone", "STING + T cells", "STING + macrophages")))]
or6
p6b = ggplot(or6, aes(x = or, y = model, colour = term)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(0.5), size = 0.4) +
  scale_x_log10() + scale_colour_manual(values = c("STING score" = "#2166AC", "T-cell score" = "#999999", "Macrophage score" = "#B2182B")) +
  theme_bw() + theme(legend.position = "bottom") +
  labs(x = "OR for metastasis per SD (log scale)", y = NULL, colour = NULL, title = "GSE21257: metastasis (34 vs 19)")
print(p6b)
ggsave("results/figures/Fig6B_STING_adjusted.pdf", p6b, width = 6, height = 3.5)
# 355. 重建 GSE21257 表达矩阵（nuID → 基因名），并保存
library(lumiHumanAll.db); library(survival)
mx  = readLines(gzfile("data_raw/gse21257/GSE21257_series_matrix.txt.gz"))
b   = grep("^!series_matrix_table_begin", mx); e = grep("^!series_matrix_table_end", mx)
ex  = fread(text = mx[(b + 1):(e - 1)])
sym = AnnotationDbi::mapIds(lumiHumanAll.db, keys = ex$ID_REF, column = "SYMBOL",
                            keytype = "PROBEID", multiVals = "first")
g21 = fread("data_proc/gse21257_clin_scores.csv")
ann = data.table(ID_REF = names(sym), gene = unname(sym))[!is.na(gene)]
ex2 = merge(ann, ex, by = "ID_REF")
ex2[, avg := rowMeans(.SD), .SDcols = g21$sample]
ex2 = ex2[order(-avg)][!duplicated(gene)]
mat21 = as.matrix(ex2[, g21$sample, with = FALSE]); rownames(mat21) = ex2$gene
saveRDS(mat21, "data_proc/gse21257_expr.rds")
dim(mat21)

# 356. 重新抽 1000 次随机基因集，保存零分布（STING–转移、T3 糖基化–生存）
z_all = t(scale(t(mat21))); z_all = z_all[complete.cases(z_all), ]
yes = g21$meta_any == "Yes"; tm = g21$os_mon; ev = g21$OS
p_meta = function(s) wilcox.test(s[yes], s[!yes])$p.value
p_cox  = function(s) summary(coxph(Surv(tm, ev) ~ scale(s)))$coefficients[5]
sting_g = intersect(c("CGAS", "STING1", "TBK1", "IRF3", "IFNB1", "CXCL10", "ISG15"), rownames(z_all))
glyc_g  = intersect(c("B4GALT1", "GALNT1", "GALNT2", "XYLT1", "B3GAT3", "MGAT1"), rownames(z_all))
set.seed(123); B = 1000
null_meta = replicate(B, p_meta(colMeans(z_all[sample(nrow(z_all), length(sting_g)), ])))
null_cox  = replicate(B, p_cox( colMeans(z_all[sample(nrow(z_all), length(glyc_g)),  ])))
obs = c(STING_meta = p_meta(colMeans(z_all[sting_g, ])), Glyco_cox = p_cox(colMeans(z_all[glyc_g, ])))
emp = c(STING_meta = mean(null_meta <= obs[1]), Glyco_cox = mean(null_cox <= obs[2]))
rbind(observed_P = signif(obs, 3), empirical_P = emp)
saveRDS(list(null_meta = null_meta, null_cox = null_cox, obs = obs, emp = emp), "results/random_null_gse21257.rds")

# 357. Fig 6C：随机基因集的 -log10(P) 分布，红线是真实基因集
nd = rbind(data.table(test = "STING score vs metastasis", lp = -log10(null_meta)),
           data.table(test = "Mn glycosylation vs survival", lp = -log10(null_cox)))
od = data.table(test = c("STING score vs metastasis", "Mn glycosylation vs survival"),
                lp = -log10(obs), lab = sprintf("observed P = %.3g\nempirical P = %.3f", obs, emp))
p6c = ggplot(nd, aes(lp)) + geom_histogram(bins = 40, fill = "grey75", colour = "white") +
  geom_vline(data = od, aes(xintercept = lp), colour = "#B2182B", linewidth = 1) +
  geom_text(data = od, aes(x = lp, y = Inf, label = lab), hjust = -0.05, vjust = 1.3, size = 3, colour = "#B2182B") +
  facet_wrap(~ test, scales = "free") + theme_bw() +
  labs(x = "-log10(P) of 1,000 random gene sets of equal size", y = "Count")
print(p6c)
ggsave("results/figures/Fig6C_random_null.pdf", p6c, width = 9, height = 3.5)
2
# ===== 检查示意图文件 =====
# 360. 列出示意图文件 ==============================
list.files("results/figures", pattern = "Fig1A|Fig7C")
# ===== 第16课：拼合组图 =====

# 361. 安装并加载拼图需要的包 ==============================
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
for (p in c("magick", "pdftools", "patchwork")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
library(magick)
library(pdftools)
library(patchwork)
library(ggplot2)
library(grid)

# 362. 盘点所有子图 PDF 及其尺寸 ==============================
pdf_files = c(list.files("results",         pattern = "\\.pdf$", full.names = TRUE),
              list.files("results/scRNA",   pattern = "\\.pdf$", full.names = TRUE),
              list.files("results/drug",    pattern = "\\.pdf$", full.names = TRUE),
              list.files("results/figures", pattern = "\\.pdf$", full.names = TRUE))
pdf_info  = data.frame(
  file    = sub("results/", "", pdf_files),
  pages   = sapply(pdf_files, function(f) pdf_info(f)$pages),
  w_inch  = sapply(pdf_files, function(f) round(pdf_pagesize(f)$width[1]  / 72, 1)),
  h_inch  = sapply(pdf_files, function(f) round(pdf_pagesize(f)$height[1] / 72, 1)),
  size_MB = round(file.size(pdf_files) / 1024^2, 1),
  row.names = NULL
)
print(pdf_info)

# 363. 定义读图函数：PDF → 300 dpi 图片块 ==============================
read_panel = function(f, page = 1, dpi = 300) {
  img  = image_read_pdf(f, pages = page, density = dpi)
  img  = image_trim(img)                                   # 去掉四周空白
  img  = image_border(img, "white", "20x20")               # 留一点边距
  wrap_elements(full = rasterGrob(as.raster(img), interpolate = TRUE))
}

# 364. 拼合 Fig 1（A 研究设计 + B MR 森林图） ==============================
fig_dir = "results/figures"
p1a     = read_panel(file.path(fig_dir, "Fig1A_study_design.pdf"))
p1b     = read_panel(file.path(fig_dir, "Fig1B_MR_forest.pdf"))
fig1    = (p1a | p1b) +
  plot_layout(widths = c(1.2, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 14, face = "bold"))
print(fig1)

# 365. 导出 Fig 1（TIFF 300 dpi + PDF） ==============================
dir.create(file.path(fig_dir, "final"), showWarnings = FALSE)
ggsave(file.path(fig_dir, "final", "Figure1.tiff"), fig1,
       width = 18, height = 9, units = "cm", dpi = 300, compression = "lzw")
ggsave(file.path(fig_dir, "final", "Figure1.pdf"),  fig1,
       width = 18, height = 9, units = "cm")
# 366. Fig 1 改为上下排列 ==============================
fig1 = (p1a / p1b) +
  plot_layout(heights = c(1.3, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 14, face = "bold"))

ggsave(file.path(fig_dir, "final", "Figure1.tiff"), fig1,
       width = 18, height = 22, units = "cm", dpi = 300, compression = "lzw")
ggsave(file.path(fig_dir, "final", "Figure1.pdf"),  fig1,
       width = 18, height = 22, units = "cm")
# 367. 确认 Figure1 已更新 ==============================
file.info(file.path(fig_dir, "final", "Figure1.tiff"))[, c("size", "mtime")]
# 368. 加大子图边距，避免标号压图 ==============================
read_panel = function(f, page = 1, dpi = 300) {
  img  = image_read_pdf(f, pages = page, density = dpi)
  img  = image_trim(img)
  img  = image_border(img, "white", "60x60")               # 边距由 20 改为 60
  wrap_elements(full = rasterGrob(as.raster(img), interpolate = TRUE))
}

# 369. 重新读图并导出 Fig 1 ==============================
p1a  = read_panel(file.path(fig_dir, "Fig1A_study_design.pdf"))
p1b  = read_panel(file.path(fig_dir, "Fig1B_MR_forest.pdf"))
fig1 = (p1a / p1b) +
  plot_layout(heights = c(1.3, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 14, face = "bold"))

ggsave(file.path(fig_dir, "final", "Figure1.tiff"), fig1,
       width = 18, height = 22, units = "cm", dpi = 300, compression = "lzw")
ggsave(file.path(fig_dir, "final", "Figure1.pdf"),  fig1,
       width = 18, height = 22, units = "cm")
# 370. 定义路径简写与导出函数 ==============================
res       = "results"
rp        = function(f, page = 1) read_panel(file.path(res, f), page)
tag_theme = theme(plot.tag = element_text(size = 14, face = "bold"))
save_fig  = function(p, name, h) {
  ggsave(file.path(fig_dir, "final", paste0(name, ".tiff")), p,
         width = 18, height = h, units = "cm", dpi = 300, compression = "lzw")
  ggsave(file.path(fig_dir, "final", paste0(name, ".pdf")),  p,
         width = 18, height = h, units = "cm")
}

# 371. 拼合 Fig 2：单细胞发现队列 ==============================
p2a  = rp("scRNA/umap_celltype.pdf")
p2b  = rp("scRNA/heatmap_mn_tiers.pdf")
p2c  = rp("scRNA/dotplot_mn_genes.pdf")
p2d  = rp("scRNA/pseudobulk_tiers_boxplot.pdf")
fig2 = ((p2a | p2b) / p2c / p2d) +
  plot_layout(heights = c(1, 0.8, 1.1)) +
  plot_annotation(tag_levels = "A") & tag_theme
save_fig(fig2, "Figure2", 24)

# 372. 拼合 Fig 3：单细胞验证 ==============================
p3a  = rp("scRNA/val_umap_celltype.pdf")
p3b  = rp("figures/Fig3C_concordance.pdf")
fig3 = (p3a | p3b) +
  plot_layout(widths = c(1, 1.6)) +
  plot_annotation(tag_levels = "A") & tag_theme
save_fig(fig3, "Figure3", 9)

# 373. 拼合 Fig 4：细胞通讯 ==============================
p4a  = rp("scRNA/cellchat_circle.pdf", page = 1)
p4b  = rp("scRNA/cellchat_bubble_mn_to_os.pdf")
p4c  = rp("scRNA/cellchat_integrin_heatmap.pdf")
p4d  = rp("scRNA/dotplot_receptors.pdf")
fig4 = ((p4a / p4c) | p4b) / p4d +
  plot_layout(heights = c(2, 0.7)) +
  plot_annotation(tag_levels = "A") & tag_theme
save_fig(fig4, "Figure4", 24)

# 374. 拼合 Fig 5：空间转录组 ==============================
p5a  = rp("spatial_example_maps_smoothed.pdf")
p5b  = rp("spatial_rho_heatmap_pub.pdf")
p5c  = rp("spatial_rho_per_patient.pdf")
p5d  = rp("spatial_tumor_primary_vs_met.pdf")
fig5 = p5a / (p5b | p5c) / p5d +
  plot_layout(heights = c(0.8, 1, 0.9)) +
  plot_annotation(tag_levels = "A") & tag_theme
save_fig(fig5, "Figure5", 22)

# 375. 拼合 Fig 6：Bulk 队列 ==============================
p6a  = rp("figures/Fig6A_STING_cell_cor.pdf")
p6b  = rp("figures/Fig6B_STING_adjusted.pdf")
p6c  = rp("figures/Fig6C_random_null.pdf")
fig6 = ((p6a | p6b) / p6c) +
  plot_annotation(tag_levels = "A") & tag_theme
save_fig(fig6, "Figure6", 14)

# 376. 拼合 Fig 7：药物转化 ==============================
p7a  = rp("drug/target_dotplot_discovery.pdf")
p7b  = rp("figures/Fig7B_trials.pdf")
p7c  = rp("figures/Fig7C_model.pdf")
fig7 = p7a / (p7b | p7c) +
  plot_layout(heights = c(1, 1.2)) +
  plot_annotation(tag_levels = "A") & tag_theme
save_fig(fig7, "Figure7", 18)
# 377. Fig 7 改版：A、B 在上，模型 C 独占下方一整行 ==============================
fig7 = (p7a / p7b / p7c) +
  plot_layout(heights = c(1, 0.55, 1.3)) +
  plot_annotation(tag_levels = "A") & tag_theme
save_fig(fig7, "Figure7", 24)
# ===== 第17课：按 JTM 规格统一作图 =====

# 378. 定义 JTM 尺寸常量与导出函数 ==============================
W_FULL    = 176                                          # 整版宽度 mm
W_HALF    = 85                                           # 半版宽度 mm
final_dir = "results/figures/final"
save_jtm  = function(p, name, h, w = W_FULL) {
  ggsave(file.path(final_dir, paste0(name, ".tiff")), p, width = w, height = h,
         units = "mm", dpi = 300, compression = "lzw")
  ggsave(file.path(final_dir, paste0(name, ".pdf")),  p, width = w, height = h,
         units = "mm", device = cairo_pdf)                 # cairo_pdf 会嵌入字体
}

# 379. 定义全文统一主题（最小 8 pt） ==============================
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

# 380. 定义细胞类型标签（去下划线） ==============================
ct_lab = c(OS_Osteoblastic   = "Osteoblastic OS",  OS_Chondroblastic = "Chondroblastic OS",
           Fibroblast_MSC    = "Fibroblast/MSC",   Myeloid           = "Myeloid",
           Osteoclast        = "Osteoclast",       Endothelial       = "Endothelial",
           Pericyte          = "Pericyte",         T_NK              = "T/NK",
           B_Plasma          = "B/plasma",         Mast              = "Mast",
           Lung_epithelial   = "Lung epithelial",  Myocyte           = "Myocyte")
tier_lab = c(T1_Transport = "T1 transport",  T2_Enzyme = "T2 Mn enzymes",
             T3_Glyco     = "T3 glycosylation", T4_STING = "T4 cGAS–STING")

# 381. 定义细胞类型统一配色（色盲友好） ==============================
ct_col = c("Osteoblastic OS"   = "#3B75AF", "Chondroblastic OS" = "#8FB8DE",
           "Fibroblast/MSC"    = "#E69F00", "Myeloid"           = "#009E73",
           "Osteoclast"        = "#CC79A7", "Endothelial"       = "#D55E00",
           "Pericyte"          = "#B8860B", "T/NK"              = "#7B3294",
           "B/plasma"          = "#C2A5CF", "Mast"              = "#999999",
           "Lung epithelial"   = "#66C2A5", "Myocyte"           = "#A6761D")

# 382. 测试主题和字体 ==============================
p_test = ggplot(data.frame(x = names(ct_col), y = 1), aes(x, y, fill = x)) +
  geom_col() + scale_fill_manual(values = ct_col) +
  labs(title = "Test: Arial 8 pt", x = NULL, y = "Value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
save_jtm(p_test, "test_theme", h = 70, w = W_HALF)
# 383. 调整 Myocyte 配色，与 Pericyte 区分 ==============================
ct_col["Myocyte"] = "#4D4D4D"
# 384. 读入单细胞对象并查看元数据 ==============================
library(Seurat)
library(data.table)
library(ggrepel)
sc   = readRDS("data_proc/sc_tme_scored.rds")
meta = sc@meta.data
print(colnames(meta))
print(table(meta$celltype))
print(Reductions(sc))

# 385. 重画 Fig 2A：UMAP（统一标签与配色） ==============================
sc$ct       = factor(unname(ct_lab[as.character(sc$celltype)]), levels = names(ct_col))
print(sum(is.na(sc$ct)))                                   # 应为 0                                 # 应为 0
umap        = as.data.table(Embeddings(sc, "umap"))
setnames(umap, c("UMAP1", "UMAP2"))
umap$ct     = sc$ct
lab_pos     = umap[, .(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2)), by = ct]
p2a = ggplot(umap, aes(UMAP1, UMAP2, colour = ct)) +
  ggrastr::rasterise(geom_point(size = 0.05, stroke = 0), dpi = 300) +
  geom_text_repel(data = lab_pos, aes(label = ct), colour = "black",
                  size = 8 / .pt, min.segment.length = 0, seed = 1) +
  scale_colour_manual(values = ct_col) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme(legend.position = "none")
save_jtm(p2a, "preview_Fig2A", h = 85, w = W_HALF)

# 386. 重画 Fig 2C：四层锰基因点图 ==============================
mn_genes = list(
  "T1 transport"     = c("SLC39A8", "SLC39A14", "SLC30A10", "SLC11A2", "ATP2C1", "ATP13A1", "TMEM165", "TFRC"),
  "T2 Mn enzymes"    = c("SOD2", "ARG2", "GLUL", "PC", "PPM1A", "PPM1B"),
  "T3 glycosylation" = c("B4GALT1", "GALNT1", "GALNT2", "XYLT1", "B3GAT3", "MGAT1"),
  "T4 cGAS–STING"    = c("CGAS", "TBK1", "IRF3", "IFNB1", "CXCL10", "ISG15"))
p2c = DotPlot(sc, features = mn_genes, group.by = "ct",
              cols = c("grey85", "#3B75AF"), dot.scale = 3.5) +
  theme_jtm() +
  theme(axis.text.x     = element_text(angle = 45, hjust = 1, face = "italic"),
        legend.position = "right") +
  labs(x = NULL, y = NULL)
save_jtm(p2c, "preview_Fig2C", h = 75)
# 387. 微调 Fig 2A：去掉标签指示线 ==============================
p2a = ggplot(umap, aes(UMAP1, UMAP2, colour = ct)) +
  ggrastr::rasterise(geom_point(size = 0.05, stroke = 0), dpi = 300) +
  geom_text_repel(data = lab_pos, aes(label = ct), colour = "black",
                  size = 8 / .pt, min.segment.length = Inf, seed = 1) +
  scale_colour_manual(values = ct_col) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme(legend.position = "none")
save_jtm(p2a, "preview_Fig2A", h = 85, w = W_HALF)

# 388. 微调 Fig 2C：图例移到底部并改标题 ==============================
p2c = p2c +
  guides(colour = guide_colourbar(title = "Scaled mean expression", title.position = "top",
                                  barwidth = unit(25, "mm"), barheight = unit(2.5, "mm")),
         size   = guide_legend(title = "% expressing", title.position = "top")) +
  theme(legend.position = "bottom",
        legend.box      = "horizontal",
        legend.margin   = margin(0, 0, 0, 0))
save_jtm(p2c, "preview_Fig2C", h = 80)
# 389. 整理评分数据（统一标签与顺序） ==============================
tiers    = c("Tier1", "Tier2", "Tier3", "Tier4")
tier_nm  = c(Tier1 = "T1 transport", Tier2 = "T2 Mn enzymes",
             Tier3 = "T3 glycosylation", Tier4 = "T4 cGAS–STING")
ct_order = c("Myeloid", "Endothelial", "Pericyte", "Osteoclast", "Fibroblast/MSC",
             "T/NK", "Chondroblastic OS", "Osteoblastic OS")
les_lab  = c(Primary = "Primary", Primary_or_Recurrent = "Primary or recurrent",
             Recurrent = "Recurrent", Lung_metastasis = "Lung metastasis")
les_col  = c("Primary" = "#0072B2", "Primary or recurrent" = "#56B4E9",
             "Recurrent" = "#E69F00", "Lung metastasis" = "#CC79A7")
dt       = as.data.table(sc@meta.data)[, c("patient", "lesion", "ct", tiers), with = FALSE]
long     = melt(dt, measure.vars = tiers, variable.name = "tier", value.name = "score")
long[, tier := factor(tier_nm[as.character(tier)], levels = tier_nm)]
long[, ct   := factor(as.character(ct), levels = rev(ct_order))]

# 390. 重画 Fig 2B：各细胞类型平均评分热图 ==============================
hm = long[, .(m = mean(score)), by = .(ct, tier)]
hm[, z := (m - mean(m)) / sd(m), by = tier]
p2b = ggplot(hm, aes(tier, ct, fill = z)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", z)), size = 8 / .pt) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, name = "z-score\n(within tier)") +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line   = element_blank(),
        axis.ticks  = element_blank())
save_jtm(p2b, "preview_Fig2B", h = 85, w = W_HALF)

# 391. 重画 Fig 2D：伪批量（每位病人）评分箱线图 ==============================
pb = long[, .(score = mean(score), n = .N), by = .(patient, lesion, ct, tier)][n >= 20]
pb[, lesion := factor(les_lab[as.character(lesion)], levels = les_lab)]
print(pb[tier == "T1 transport", .N, by = ct][order(ct)])   # 每种细胞纳入的病人数
p2d = ggplot(pb, aes(score, ct)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.3, width = 0.6, fill = NA) +
  geom_point(aes(colour = lesion), size = 0.8,
             position = position_jitter(height = 0.15, seed = 1)) +
  facet_wrap(~ tier, nrow = 1, scales = "free_x") +
  scale_colour_manual(values = les_col, name = "Lesion") +
  labs(x = "Mean module score per patient", y = NULL) +
  theme(legend.position    = "bottom",
        panel.grid.major.x = element_line(linewidth = 0.2, colour = "grey90"))
save_jtm(p2d, "preview_Fig2D", h = 70)

# 392. 拼合并导出 Fig 2 ==============================
fig2 = (p2a | p2b) / p2c / p2d +
  plot_layout(heights = c(1.1, 1, 1)) +
  plot_annotation(tag_levels = "A")
save_jtm(fig2, "Figure2", h = 220)
# 393. 修正 Fig 2：A 取消对齐，D 减少刻度 ==============================
p2d = p2d + scale_x_continuous(breaks = scales::breaks_pretty(n = 3))
fig2 = (free(p2a) | p2b) / p2c / p2d +
  plot_layout(heights = c(1.1, 1, 1)) +
  plot_annotation(tag_levels = "A")
save_jtm(fig2, "Figure2", h = 220)
# 394. 保存 Fig 2 子图对象 ==============================
saveRDS(list(a = p2a, b = p2b, c = p2c, d = p2d),
        "results/figures/final/fig2_panels.rds")

# 395. 查看 Fig 3 所需数据 ==============================
sc_val = readRDS("data_proc/sc_val_annotated.rds")
print(table(sc_val$celltype))
print(Reductions(sc_val))
conc = fread("results/scRNA/discovery_validation_combined.csv")
print(colnames(conc))
print(conc)
# 394. 保存 Fig 2 子图对象（修正版） ==============================
f_rds = "results/figures/final/fig2_panels.rds"
if (file.exists(f_rds)) file.remove(f_rds)
for (nm in c("p2a", "p2b", "p2c", "p2d")) {
  p          = get(nm)
  p$plot_env = globalenv()                                 # 断开与 sc 的引用
  assign(nm, p)
}
saveRDS(list(a = p2a, b = p2b, c = p2c, d = p2d), f_rds)
print(round(file.size(f_rds) / 1024^2, 1))                 # 文件大小 MB
# 394. 保存 Fig 2 绘图数据（替代保存图对象） ==============================
f_rds = "results/figures/final/fig2_data.rds"
if (file.exists("results/figures/final/fig2_panels.rds"))
  file.remove("results/figures/final/fig2_panels.rds")
saveRDS(list(umap = umap, lab_pos = lab_pos, hm = hm, pb = pb, dot = p2c$data), f_rds)
print(round(file.size(f_rds) / 1024^2, 1))                 # 文件大小 MB
# 396. 定义 UMAP 绘图函数（Fig 2A / 3A 共用） ==============================
plot_umap = function(obj) {
  d   = as.data.table(Embeddings(obj, "umap"))
  setnames(d, c("UMAP1", "UMAP2"))
  d$ct = obj$ct
  lab = d[, .(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2)), by = ct]
  ggplot(d, aes(UMAP1, UMAP2, colour = ct)) +
    ggrastr::rasterise(geom_point(size = 0.05, stroke = 0), dpi = 300) +
    geom_text_repel(data = lab, aes(label = ct), colour = "black",
                    size = 8 / .pt, min.segment.length = Inf, seed = 1) +
    scale_colour_manual(values = ct_col) +
    labs(x = "UMAP 1", y = "UMAP 2") +
    theme(legend.position = "none")
}

# 397. 重画 Fig 3A：验证队列 UMAP ==============================
sc_val$ct = factor(unname(ct_lab[as.character(sc_val$celltype)]), levels = names(ct_col))
print(sum(is.na(sc_val$ct)))                               # 应为 0
p3a = plot_umap(sc_val)

# 398. 整理一致性数据 ==============================
tier_short = c(T1_Transport = "Mn transport", T2_Enzyme = "Mn enzymes",
               T3_Glyco     = "Mn glycosylation", T4_STING = "cGAS–STING")
fmt_p      = function(p) ifelse(p < 0.001, sprintf("%.1e", p), sprintf("%.3f", p))
conc[, lab := paste0(ct_lab[cell], " — ", tier_short[tier])]
conc       = conc[order(p_sign_BH, -k_all / n_all)]
conc[, lab := factor(lab, levels = rev(lab))]
bar        = rbind(conc[, .(lab, cohort = "Discovery (n = 11)",  k = k_disc, n = n_disc)],
                   conc[, .(lab, cohort = "Validation (n = 6)",  k = k_val,  n = n_val)])
bar[, prop   := k / n]
bar[, cohort := factor(cohort, levels = c("Validation (n = 6)", "Discovery (n = 11)"))]

# 399. 重画 Fig 3B：两队列方向一致性条形图 ==============================
p3b = ggplot(bar, aes(prop, lab, fill = cohort)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  geom_text(aes(label = paste0(k, "/", n)), position = position_dodge(width = 0.75),
            hjust = 1.15, size = 8 / .pt, colour = "white") +
  geom_vline(xintercept = 0.5, linetype = "dashed", linewidth = 0.3, colour = "grey40") +
  geom_text(data = conc, inherit.aes = FALSE, hjust = 0, size = 8 / .pt,
            aes(x = 1.04, y = lab, label = paste0(k_all, "/", n_all, "   ", fmt_p(p_sign_BH)))) +
  annotate("text", x = 1.04, y = Inf, label = "Pooled   BH P", hjust = 0, vjust = -0.6,
           size = 8 / .pt, fontface = "bold") +
  scale_x_continuous(labels = scales::percent, breaks = c(0, 0.5, 1), expand = c(0, 0)) +
  scale_fill_manual(values = c("Discovery (n = 11)" = "grey30", "Validation (n = 6)" = "grey65"),
                    breaks = c("Discovery (n = 11)", "Validation (n = 6)"), name = NULL) +
  coord_cartesian(xlim = c(0, 1), clip = "off") +
  labs(x = "Patients with concordant direction\n(vs osteoblastic OS cells)", y = NULL) +
  theme(legend.position = "bottom",
        plot.margin     = margin(15, 75, 5, 5))

# 400. 拼合并导出 Fig 3 ==============================
fig3 = (free(p3a) | p3b) +
  plot_layout(widths = c(0.8, 1.2)) +
  plot_annotation(tag_levels = "A")
save_jtm(fig3, "Figure3", h = 95)
saveRDS(list(conc = conc, bar = bar), "results/figures/final/fig3_data.rds")
# 401. 改用 ragg 导出 TIFF ==============================
if (!requireNamespace("ragg", quietly = TRUE)) install.packages("ragg")
save_jtm = function(p, name, h, w = W_FULL) {
  ggsave(file.path(final_dir, paste0(name, ".tiff")), p, width = w, height = h,
         units = "mm", dpi = 300, device = ragg::agg_tiff, compression = "lzw")
  ggsave(file.path(final_dir, paste0(name, ".pdf")),  p, width = w, height = h,
         units = "mm", device = cairo_pdf)
}

# 402. 重新导出 Fig 2 和 Fig 3 并检查时间 ==============================
save_jtm(fig2, "Figure2", h = 220)
save_jtm(fig3, "Figure3", h = 95)
file.info(file.path(final_dir, c("Figure2.tiff", "Figure3.tiff")))[, c("size", "mtime")]
# 403. UMAP 标签加白色描边 ==============================
plot_umap = function(obj) {
  d   = as.data.table(Embeddings(obj, "umap"))
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

# 404. 重新导出 Fig 3（加宽 A） ==============================
p3a  = plot_umap(sc_val)
fig3 = (free(p3a) | p3b) +
  plot_layout(widths = c(1.5, 1)) +
  plot_annotation(tag_levels = "A")
save_jtm(fig3, "Figure3", h = 95)

# 405. Fig 2A 同步使用新函数 ==============================
p2a  = plot_umap(sc)
fig2 = (free(p2a) | p2b) / p2c / p2d +
  plot_layout(heights = c(1.1, 1, 1)) +
  plot_annotation(tag_levels = "A")
save_jtm(fig2, "Figure2", h = 220)
ggsave(file.path(final_dir, "Figure3.tiff"), fig3, width = 176, height = 95,
       units = "mm", dpi = 300, device = ragg::agg_tiff, compression = "lzw")
file.info(file.path(final_dir, c("Figure2.tiff", "Figure3.tiff", "Figure3.pdf")))[, c("size", "mtime")]
# 407. 加载全局主题与函数 ==============================
source("scripts/15_figure_theme.R")               # 路径按你保存的位置修改
library(CellChat)
library(data.table)
library(patchwork)

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
# ===== 15_figure_theme.R：JTM 统一作图规格 =====

# 1. 加载作图相关包 ==============================
library(ggplot2)
library(patchwork)
library(data.table)
library(ggrepel)

# 2. 尺寸常量与输出目录 ==============================
W_FULL    = 176                                          # 整版宽度 mm
W_HALF    = 85                                           # 半版宽度 mm
final_dir = "results/figures/final"
dir.create(final_dir, showWarnings = FALSE, recursive = TRUE)

# 3. 导出函数（TIFF 300 dpi + 嵌入字体 PDF） ==============================
save_jtm = function(p, name, h, w = W_FULL) {
  ggsave(file.path(final_dir, paste0(name, ".tiff")), p, width = w, height = h,
         units = "mm", dpi = 300, device = ragg::agg_tiff, compression = "lzw")
  ggsave(file.path(final_dir, paste0(name, ".pdf")),  p, width = w, height = h,
         units = "mm", device = cairo_pdf)
}

# 4. 全文统一主题（最小 8 pt） ==============================
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

# 5. 细胞类型与基因层级标签（去下划线） ==============================
ct_lab = c(OS_Osteoblastic   = "Osteoblastic OS",  OS_Chondroblastic = "Chondroblastic OS",
           Fibroblast_MSC    = "Fibroblast/MSC",   Myeloid           = "Myeloid",
           Osteoclast        = "Osteoclast",       Endothelial       = "Endothelial",
           Pericyte          = "Pericyte",         T_NK              = "T/NK",
           B_Plasma          = "B/plasma",         Mast              = "Mast",
           Lung_epithelial   = "Lung epithelial",  Myocyte           = "Myocyte")
tier_lab = c(T1_Transport = "T1 transport",     T2_Enzyme = "T2 Mn enzymes",
             T3_Glyco     = "T3 glycosylation", T4_STING  = "T4 cGAS–STING")

# 6. 细胞类型统一配色（色盲友好） ==============================
ct_col = c("Osteoblastic OS"   = "#3B75AF", "Chondroblastic OS" = "#8FB8DE",
           "Fibroblast/MSC"    = "#E69F00", "Myeloid"           = "#009E73",
           "Osteoclast"        = "#CC79A7", "Endothelial"       = "#D55E00",
           "Pericyte"          = "#B8860B", "T/NK"              = "#7B3294",
           "B/plasma"          = "#C2A5CF", "Mast"              = "#999999",
           "Lung epithelial"   = "#66C2A5", "Myocyte"           = "#4D4D4D")

# 7. P 值格式化函数 ==============================
fmt_p = function(p) ifelse(p < 0.001, sprintf("%.1e", p), sprintf("%.3f", p))

# 8. UMAP 绘图函数（带白色描边标签） ==============================
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
  # 412. 检查 data.table 是否完好 ==============================
  library(data.table)
  packageVersion("data.table")
  as.data.table(data.frame(a = 1:3))
}
# 413. 检查主题脚本中是否含安装命令 ==============================
grep("install|update", readLines("scripts/15_figure_theme.R"), value = TRUE)
length(readLines("scripts/15_figure_theme.R"))
# 414. 检查 17:00 后是否有数据或结果文件被改动 ==============================
fs = list.files(c("data_proc", "results"), recursive = TRUE, full.names = TRUE)
fi = file.info(fs)
fi[fi$mtime > as.POSIXct("2026-11-23 17:00:00"), c("size", "mtime")]
