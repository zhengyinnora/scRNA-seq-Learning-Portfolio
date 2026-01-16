# Load essential libraries for single-cell analysis pipeline
library(Seurat)    # Core framework for scRNA-seq analysis and visualization
library(dplyr)     # Data manipulation and dataframe filtering
library(patchwork) # Multi-plot composition and layout arrangement
library(Matrix)    # Efficient handling of sparse matrices for memory optimization

# ==============================================================================
# Step 1: Verify Working Directory
# [EN] Check if the current working directory corresponds to the project root.
#      This ensures relative paths work correctly across different environments.
# [CN] 检查当前工作目录是否已定位至项目根目录。
#      此步骤用于确保相对路径在不同运行环境下均可正确索引。
# ==============================================================================
print(getwd())

# ==============================================================================
# Step 2: Define Directory Paths
# [EN] Set up input and output paths relative to the project root.
#      'raw_data_dir' stores the original matrix files.
#      'output_dir' is designated for saving processed results.
# [CN] 定义相对文件路径。
#      'raw_data_dir' 用于存放原始矩阵文件；
#      'output_dir' 指定为处理后结果的保存路径。
# ==============================================================================

# 定义原始数据路径 (Input)
raw_data_dir <- "lianxi/00_data_raw"

# 定义输出路径 (Output)
output_dir <- "lianxi/01_data_processed"

# ==============================================================================
# Validation: Check file visibility
# [EN] Verify if the raw data files are visible in the specified directory.
# [CN] 验证指定路径下的原始数据文件是否可见。
# ==============================================================================

list.files(raw_data_dir)

# ==============================================================================
# Step 3: Load Sparse Matrix Data
# [EN] Load the 10X Genomics data by mapping the matrix, barcodes, and features files.
#      'ReadMtx' is used here because the file names have custom prefixes (GSE...).
# [CN] 读取 10X Genomics 稀疏矩阵数据。
#      通过映射矩阵(mtx)、条形码(cells)和基因特征(features)文件来组装数据。
#      此处使用 'ReadMtx' 函数以处理带有自定义前缀的文件名。
# ==============================================================================

counts <- ReadMtx(
  mtx      = file.path(raw_data_dir, "GSE141259_WholeLung_rawcounts.mtx.gz"),
  cells    = file.path(raw_data_dir, "GSE141259_WholeLung_barcodes.txt.gz"),
  features = file.path(raw_data_dir, "GSE141259_WholeLung_genes.txt.gz"),
  feature.column = 1  # [CN] 指定基因名在文件的第1列 / [EN] Gene names are in column 1
)

# ==============================================================================
# Validation: Check Data Dimensions
# [EN] Inspect the dimensions of the loaded matrix (Genes x Cells).
# [CN] 检查加载后的矩阵维度（基因数 x 细胞数）。
# ==============================================================================

dim(counts)

# ==============================================================================
# Step 4: Initialize Seurat Object
# [EN] Create the Seurat object. Note: Seurat auto-parses cell barcodes (e.g., 'muc4169_...')
#      to populate the 'orig.ident' column, often overriding the 'project' argument.
# [CN] 初始化 Seurat 对象。注意：Seurat 会自动解析细胞条形码（如 'muc4169_...'）
#      来填充 'orig.ident' 列，这通常会覆盖我们设置的 'project' 参数。
# ==============================================================================

lung_obj <- CreateSeuratObject(
  counts = counts,
  project = "lianxixi", 
  min.cells = 3,
  min.features = 200
)

# ==============================================================================
# Step 4.1: Standardize Metadata Identity
# [EN] Standardization of 'orig.ident'.
#      Seurat automatically infers sample identity from cell barcodes (e.g., 'muc4169_...').
#      Here, we unify the identity of all cells to the project name "lianxixi" 
#      to ensure consistent grouping for downstream visualization.
# [CN] 标准化元数据身份标签。
#      Seurat 默认会从细胞条形码中推断样本身份（如 'muc4169_...'）。
#      此处我们将所有细胞的身份标签统一修正为项目名称 "lianxixi"，
#      以确保后续可视化分析时的分组一致性。
# ==============================================================================

lung_obj$orig.ident <- "lianxixi"

# [Check] 再次检查，现在必须是 "lianxixi" 了
head(lung_obj@meta.data)

# ==============================================================================
# Step 5: Load External Metadata
# [EN] Load the metadata CSV file which contains cell type annotations provided by the authors.
#      Note: We set 'row.names = 1' to use the first column (cell barcodes) as row identifiers.
# [CN] 读取外部元数据文件（包含作者提供的细胞类型注释）。
#      注意：设置 'row.names = 1' 以将第一列（细胞条形码）作为行索引。
# ==============================================================================

# 读取 CSV 文件
metadata_csv <- read.csv(file.path(raw_data_dir, "GSE141259_WholeLung_cellinfo.csv.gz"), row.names = 1)

# ==============================================================================
# Validation: Inspect Barcode Formats
# [EN] Critical Check: Compare the row names of the loaded metadata with the column names of the Seurat object.
#      They MUST match exactly for the metadata to be added correctly.
# [CN] 关键检查：对比新读取的元数据行名与 Seurat 对象的列名。
#      两者必须完全一致（包括标点符号），才能正确添加元数据。
# ==============================================================================

print("--- Seurat Object Cell Names (My Data) ---")
head(colnames(lung_obj))

print("--- Author's Metadata Cell Names (Answer Key) ---")
head(rownames(metadata_csv))

# ==============================================================================
# Step 6: Integrate Metadata
# [EN] Add the loaded metadata table to the Seurat object.
#      Since the cell names matched perfectly, Seurat will automatically map
#      each row of the metadata to the corresponding cell in the object.
# [CN] 整合元数据。
#      将读取的元数据表添加到 Seurat 对象中。
#      由于细胞名称完全匹配，Seurat 会自动将元数据里的每一行信息映射到对应的细胞上。
# ==============================================================================

lung_obj <- AddMetaData(lung_obj, metadata = metadata_csv)

# ==============================================================================
# Validation: Inspect New Columns
# [EN] Check the metadata table again. You should see new columns (e.g., cell_type) added to the right.
# [CN] 再次检查元数据表。你应该能看到右侧新增了包含作者注释信息的列（如 cell_type）。
# ==============================================================================

head(lung_obj@meta.data)

# ==============================================================================
# Step 7: Standard Processing Pipeline
# [EN] Run the standard Seurat workflow to visualize the data.
#      1. NormalizeData: Standardize gene expression across cells.
#      2. FindVariableFeatures: Identify highly variable genes (HVGs) for analysis.
#      3. ScaleData: Shift/scale expression so that highly-expressed genes don't dominate.
#      4. RunPCA: Dimensionality reduction (compress data to ~50 dims).
#      5. RunUMAP: Non-linear embedding to visualize clusters in 2D space.
# [CN] 运行标准 Seurat 处理流程以进行可视化。
#      1. NormalizeData: 对细胞间的基因表达量进行标准化。
#      2. FindVariableFeatures: 寻找高变基因（HVG）用于后续分析。
#      3. ScaleData: 对数据进行中心化和缩放，防止高表达基因主导分析结果。
#      4. RunPCA: 降维分析（将数据压缩至主成分）。
#      5. RunUMAP: 非线性降维，将细胞群投射到二维平面进行可视化。
# ==============================================================================

# 1. 标准化
lung_obj <- NormalizeData(lung_obj, normalization.method = "LogNormalize", scale.factor = 10000)

# 2. 找高变基因 (找前 2000 个差异最大的基因)
lung_obj <- FindVariableFeatures(lung_obj, selection.method = "vst", nfeatures = 2000)

# 3. 归一化 (这一步可能稍慢，耐心等待)
lung_obj <- ScaleData(lung_obj)

# 4. PCA 降维
lung_obj <- RunPCA(lung_obj, features = VariableFeatures(object = lung_obj))

# 5. UMAP 可视化计算 (这一步可能需要几分钟)
lung_obj <- RunUMAP(lung_obj, dims = 1:30)

# ==============================================================================
# Step 8 (Fix NA): Handle Missing Labels
# [EN] Diagnosis: The error "differing number of rows" occurs because 29 cells 
#      have 'NA' (missing values) in the 'cell.type' column.
#      Fix: We assign the label "Unknown" to these cells to prevent plotting errors.
# [CN] 诊断：报错显示行数不一致，是因为有 29 个细胞的 'cell.type' 列是 NA（缺失值）。
#      修复：我们将这些缺失值的细胞标记为 "Unknown"，以确保绘图函数正常运行。
# ==============================================================================

# 1. 把这一列先转成字符格式（防止因为是 Factor 因子而报错）
lung_obj$cell.type <- as.character(lung_obj$cell.type)

# 2. 找到那些是 NA 的细胞，给它们起名叫 "Unknown"
lung_obj$cell.type[is.na(lung_obj$cell.type)] <- "Unknown"

# ==============================================================================
# Final Attempt: Visualization
# [EN] Now that all cells have a label, the plot should render correctly.
# [CN] 现在所有细胞都有了标签，UMAP 图应该可以完美呈现了。
# ==============================================================================

DimPlot(lung_obj, reduction = "umap", group.by = "cell.type", label = TRUE)

# ==============================================================================
# Step 8 (Refined): The "Open & Generous" Plot
# [CN] 优化绘图：
#      1. label = TRUE: 在图上直接显示名字。
#      2. repel = TRUE: 让文字自动避让，防止重叠（关键参数！）。
#      3. NoLegend(): 关掉右边那个巨大的图例列表，让图像铺满屏幕。
# ==============================================================================

DimPlot(lung_obj, reduction = "umap", group.by = "cell.type", label = TRUE, repel = TRUE) + NoLegend()

# ==============================================================================
# Step 9 (Fixed): Save to Specific Folder
# [CN] 修正：将高清大图保存到指定的 "04_output_plots" 文件夹。
# [EN] Fix: Save the high-resolution plot to the specific "04_output_plots" directory.
# ==============================================================================

# 1. 定义图片保存的路径 (指向那个 04 文件夹)
plots_dir <- "lianxi/04_output_plots"

# 2. 再次把图赋值给变量 p1 (以防万一)
p1 <- DimPlot(lung_obj, reduction = "umap", group.by = "cell.type", label = TRUE, repel = TRUE) + NoLegend()

# 3. 使用 file.path 拼接完整的路径进行保存
# 最终路径会是: lianxi/04_output_plots/lianxixi_UMAP_Generous.png
ggsave(filename = file.path(plots_dir, "lianxixi_UMAP_Generous.png"), plot = p1, width = 12, height = 10)

# ==============================================================================
# Step 10: Save Processed Object (Checkpoint)
# [CN] 保存处理后的 Seurat 对象。
#      这是一个关键的“存档点”。将包含 UMAP 坐标和元数据的最终对象保存为 .rds 文件。
#      下次分析时，只需使用 readRDS() 读取此文件，无需重新运行之前的耗时步骤。
# [EN] Save the processed Seurat object as an .rds file.
#      This acts as a checkpoint, preserving all results (UMAP, metadata).
#      Use readRDS() to reload this object instantly in future sessions.
# ==============================================================================

# 定义保存路径
save_path <- file.path(output_dir, "lung_obj_processed_with_celltype.rds")

# 保存对象 (这可能需要几十秒，视文件大小而定)
saveRDS(lung_obj, file = save_path)

# 打印提示
message("Saved successfully to: ", save_path)

