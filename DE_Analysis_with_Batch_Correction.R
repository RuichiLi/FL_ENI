# 加载必要的包
library(readxl)   # 用于读取Excel文件
library(sva)      # 用于批次效应校正
library(limma)    # 用于差异表达分析
library(edgeR)    # 用于count数据标准化
library(openxlsx) # 用于输出Excel结果

# 设置工作目录（替换为实际路径）
setwd("H:/FL_ENI")

# 读取原始count数据（Excel的第一个sheet）
count_data <- read_excel("RNAseq.xlsx", sheet = 1)

# 将第一列转换为行名
count_matrix <- as.matrix(count_data[, -1])
rownames(count_matrix) <- count_data[[1]]

# 读取批次信息（Excel的第二个sheet）
batch_info <- read_excel("RNAseq.xlsx", sheet = 2)
# 确保样本顺序一致
batch_vector <- batch_info$Batch[match(colnames(count_matrix), batch_info$Sample)]

# 读取分组信息（Excel的第三个sheet）
group_info <- read_excel("RNAseq.xlsx", sheet = 3)
# 确保样本顺序一致
group_vector <- factor(group_info$Group[match(colnames(count_matrix), group_info$Sample)])

# 检查数据一致性
if(!all(colnames(count_matrix) == group_info$Sample[match(colnames(count_matrix), group_info$Sample)])) {
  stop("样本名在count矩阵和分组信息中不匹配！")
}

# Step 1: 数据标准化 - 使用edgeR的TMM标准化
dge <- DGEList(counts = count_matrix)
dge <- calcNormFactors(dge, method = "TMM")

# Step 2: 转换为log2-CPM值
logCPM <- cpm(dge, log = TRUE, prior.count = 3)

# Step 3: 使用ComBat进行批次效应校正
design <- model.matrix(~group_vector)
corrected_data <- ComBat(
  dat = logCPM,
  batch = batch_vector,
  mod = design,
  par.prior = TRUE,
  prior.plots = FALSE
)

# Step 4: 差异表达分析
design_diff <- model.matrix(~0 + group_vector)
colnames(design_diff) <- levels(group_vector)

fit <- lmFit(corrected_data, design_diff)

# 自动创建对比组（假设只有两个组别）
if(nlevels(group_vector) == 2) {
  groups <- levels(group_vector)
  contrast_name <- paste(groups[2], groups[1], sep = "-")
  contrast_matrix <- makeContrasts(
    contrasts = contrast_name,
    levels = design_diff
  )
} else {
  stop("当前代码仅支持两组比较，请手动设置对比组")
}

fit_contrasts <- contrasts.fit(fit, contrast_matrix)
fit_contrasts <- eBayes(fit_contrasts)

# 提取所有差异表达结果
deg_results <- topTable(
  fit_contrasts,
  coef = 1,
  number = Inf,
  adjust.method = "BH"
)

# 添加基因名列
deg_results$Gene <- rownames(deg_results)

# 保存结果到Excel
output <- createWorkbook()

# 添加DEG结果
addWorksheet(output, "DEG_Results")
writeData(output, sheet = "DEG_Results", deg_results)

# 添加校正后的表达矩阵
corrected_df <- as.data.frame(corrected_data)
corrected_df$Gene <- rownames(corrected_data)
addWorksheet(output, "Corrected_Expression")
writeData(output, sheet = "Corrected_Expression", corrected_df)

# 添加分析参数信息
params <- data.frame(
  Parameter = c("Analysis Date", "R Version", "sva Version", "limma Version", 
                "Prior Count", "Batch Correction Method", "Contrast"),
  Value = c(as.character(Sys.Date()), 
            R.version.string, 
            as.character(packageVersion("sva")),
            as.character(packageVersion("limma")),
            "3",
            "ComBat",
            contrast_name)
)
addWorksheet(output, "Analysis_Parameters")
writeData(output, sheet = "Analysis_Parameters", params)

# 保存Excel文件
saveWorkbook(output, "Analysis_Results.xlsx", overwrite = TRUE)

# 生成结果摘要
cat("差异表达分析完成！\n")
cat("总基因数:", nrow(deg_results), "\n")
cat("显著差异基因数 (adj.P.Val < 0.05):", sum(deg_results$adj.P.Val < 0.05), "\n")
cat("结果已保存至 Analysis_Results.xlsx\n")
cat("包含三个工作表:\n")
cat("1. DEG_Results: 差异表达基因结果\n")
cat("2. Corrected_Expression: 批次校正后的表达矩阵\n")
cat("3. Analysis_Parameters: 分析参数记录\n")