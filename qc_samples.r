# Carregar pacotes necessários
library(Seurat)
library(SeuratWrappers)
library(hdf5r)
library(dplyr)
library(ggplot2)
library(scDblFinder)
library(SingleCellExperiment)
library(patchwork)

setwd("C:/Users/joyce/OneDrive - FIOCRUZ/doutorado")

# -------------------------------
# 1. Defina todas as amostras
# -------------------------------

# --- lista DENV ---
arquivos_h5 <- list(
  dt1_control = "dataset1/denv_control/filtered_feature_bc_matrix.h5",
  dt1_DF = "dataset1/denv_DF/filtered_feature_bc_matrix.h5",
  dt1_DHF = "dataset1/denv_DHF/filtered_feature_bc_matrix.h5",
  dt2_DF_1  = "dataset2/denv_nat01/filtered_feature_bc_matrix.h5",
  dt2_DF_2 = "dataset2/denv_nat02/filtered_feature_bc_matrix.h5",
  dt3_primary = "dataset3/denv_primary/filtered_feature_bc_matrix.h5",
  dt3_secundary = "dataset3/denv_secundary/filtered_feature_bc_matrix.h5",
  dt4_control = "dataset4/denv_control/filtered_feature_bc_matrix.h5",
  dt4_DF = "dataset4/denv_DF/filtered_feature_bc_matrix.h5",
  dt4_DWS = "dataset4/denv_DWS/filtered_feature_bc_matrix.h5",
  dt4_SD = "dataset4/denv_SD/filtered_feature_bc_matrix.h5"
)

# --- barcode DENV ---
amostras_por_arquivo <- list(
  dt1_control = "Healthy_Control_run1",
  dt1_DF = c("DF_Day_minus_1_run1","DF_Day_minus_1_run2","DF_Day_minus_2_run1","DF_Def_run1",
            "DF_Def_run2","DF_Wk2_run1"),
  dt1_DHF = c("DHF_Day_minus_1_run1","DHF_Day_minus_1_run2","DHF_Day_minus_2_run1",
              "DHF_Def_run1","DHF_Def_run2","DHF_Wk2_run1"),
  dt2_DF_1  = c("SRR12215051","SRR12215052","SRR12215053"),
  dt2_DF_2 = c("SRR12215054","SRR12215055","SRR12215056"),
  dt3_primary = c("SRR11088622_Primary1_D1","SRR11088623_Primary1_D3","SRR11088624_Primary2_D1",
                 "SRR11088625_Primary2_D5","SRR11088626_Primary3_D1","SRR11088627_Primary3_D2"),
  dt3_secundary = c("SRR11088628_Secondary1_D0","SRR11088629_Secondary1_D1","SRR11088630_Secondary1_D7",
                 "SRR11088631_Secondary2_D1","SRR11088632_Secondary2_D1","SRR11088633_Secondary2_D5",
                 "SRR11088634_Secondary3_D0","SRR11088635_Secondary3_D1","SRR11088636_Secondary3_D5"),
  dt4_control = c("SRR22739533","SRR22739534","SRR22739543","SRR22739544","SRR22739551","SRR22739552"),
  dt4_DF = c("SRR22739530","SRR22739535","SRR22739536","SRR22739537","SRR22739539","SRR22739541",
            "SRR22739545","SRR22739547","SRR22739549","SRR22739550","SRR22739555"),
  dt4_DWS = c("SRR22739525","SRR22739527","SRR22739529","SRR22739538"),
  dt4_SD = c("SRR22739526","SRR22739528","SRR22739531","SRR22739532","SRR22739540","SRR22739542",
             "SRR22739546","SRR22739548","SRR22739553","SRR22739554")
)

# --- Lista final ---
seurat_list <- list()

# -------------------------------
# 2. Processar cada arquivo .h5
# -------------------------------

for (nome_arquivo in names(arquivos_h5)) {
  caminho <- arquivos_h5[[nome_arquivo]]
  amostras <- amostras_por_arquivo[[nome_arquivo]]
  
  # Lê o arquivo
  counts <- Read10X_h5(caminho)
  
  total_cells <- ncol(counts)
  n_amostras <- length(amostras)
  cells_por_amostra <- floor(total_cells / n_amostras)
  
  for (i in seq_along(amostras)) {
    start <- (i - 1) * cells_por_amostra + 1
    end <- if (i == n_amostras) total_cells else i * cells_por_amostra
    
    # Seleciona subconjunto de células
    barcodes <- colnames(counts)[start:end]
    sub_counts <- counts[, barcodes]
    
    # Renomear as células com prefixo da amostra
    colnames(sub_counts) <- paste0(amostras[i], "_", barcodes)
    
    # Criar objeto Seurat
    seu <- CreateSeuratObject(counts = sub_counts, project = amostras[i])
    seu$amostra <- amostras[i]
    
    seurat_list[[amostras[i]]] <- seu
  }
}
# -------------------------------
# 3. Unir todos os objetos em um único
# -------------------------------

# --- Unir todos os objetos em um só ---
combined <- merge(
  seurat_list[[1]],
  y = seurat_list[-1],
  add.cell.ids = names(seurat_list),
  project = "TodasAmostras"
)

# -------------------------------
# 4. Cálculo do %MT e criação dos violins
# -------------------------------

combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^MT-")

combined$grupo <- case_when(
  combined$amostra %in% c("Healthy_Control_run1") ~ "dt1_control",
  combined$amostra %in% c("DF_Day_minus_1_run1","DF_Day_minus_1_run2","DF_Day_minus_2_run1",
                        "DF_Def_run1","DF_Def_run2,DF_Wk2_run1") ~ "dt1_DF",
  combined$amostra %in% c("DHF_Day_minus_1_run1","DHF_Day_minus_1_run2","DHF_Day_minus_2_run1","DHF_Def_run1",
                        "DHF_Def_run2","DHF_Wk2_run1") ~ "dt1_DHF",
  combined$amostra %in% c("SRR12215051","SRR12215052","SRR12215053") ~ "dt2_DF_1",
  combined$amostra %in% c("SRR12215054","SRR12215055","SRR12215056") ~ "dt2_DF_2",
  combined$amostra %in% c("SRR11088622_Primary1_D1","SRR11088623_Primary1_D3","SRR11088624_Primary2_D1",
                        "SRR11088625_Primary2_D5","SRR11088626_Primary3_D1","SRR11088627_Primary3_D2") ~ "dt3_primary",
  combined$amostra %in% c("SRR11088628_Secondary1_D0","SRR11088629_Secondary1_D1","SRR11088630_Secondary1_D7",
                        "SRR11088631_Secondary2_D1","SRR11088632_Secondary2_D1","SRR11088633_Secondary2_D5",
                        "SRR11088634_Secondary3_D0","SRR11088635_Secondary3_D1","SRR11088636_Secondary3_D5") ~ "dt3_secundary",
  combined$amostra %in% c("SRR22739533","SRR22739534","SRR22739543","SRR22739544","SRR22739551","SRR22739552") ~ "dt4_control",
  combined$amostra %in% c("SRR22739530","SRR22739535","SRR22739536","SRR22739537","SRR22739539","SRR22739541",
                          "SRR22739545","SRR22739547","SRR22739549","SRR22739550","SRR22739555") ~ "dt4_DF",
  combined$amostra %in% c("SRR22739525","SRR22739527","SRR22739529","SRR22739538") ~ "dt4_DWS",
  combined$amostra %in% c("SRR22739526","SRR22739528","SRR22739531","SRR22739532","SRR22739540","SRR22739542",
             "SRR22739546","SRR22739548","SRR22739553","SRR22739554") ~ "dt4_SD",
  TRUE ~ "outros"
)


group_colors <- c(
  # dt1 → tons de vermelho
  "dt1_control" = "black",   # vermelho escuro
  "dt1_DF"      = "#de2d26",   # vermelho médio
  "dt1_DHF"     = "#fb6a4a",   # vermelho claro

  "dt2_DF_1"    = "#08519c",   # azul escuro
  "dt2_DF_2"    = "#3182bd",   # azul médio

  "dt3_primary"   = "#238b45", # verde escuro
  "dt3_secundary" = "#74c476", # verde claro

  "dt4_control" = "black",   # roxo escuro
  "dt4_DF"      = "#6a51a3",   # roxo médio
  "dt4_DWS"     = "#ae017e",   # rosa
  "dt4_SD"      = "#dd3497",   # rosa claro

  "outros"      = "#999999"    # cinza neutro
)


# --- Extrair metadados e valores ---
df <- FetchData(combined, vars = c("nFeature_RNA", "amostra", "grupo"))
df_qc <- FetchData(combined, vars = c("nFeature_RNA", "nCount_RNA", "amostra", "grupo"))
amostras <- c("dt1", "dt2", "dt3", "dt4")
df_nCount   <- FetchData(combined, vars = c("nCount_RNA", "amostra", "grupo"))
df_percent  <- FetchData(combined, vars = c("percent.mt", "amostra", "grupo"))

# --- cria diretório ---
output_dir <- "violin_plots_all"
dir.create(output_dir, showWarnings = FALSE)

# --- construir violin plots ---

# --- nFeature ---
vln_nFeature <- ggplot(df, aes(x = amostra, y = nFeature_RNA, fill = grupo)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_jitter(width = 0.2, size = 0.2, alpha = 0.1) +
  scale_fill_manual(values = group_colors) +
  ylim(0, 8000) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "nFeature_RNA", y = " ", x = " ")

vln_nFeature_group <- ggplot(df, aes(x = grupo, y = nFeature_RNA, fill = grupo)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_jitter(width = 0.2, size = 0.2, alpha = 0.1) +
  scale_fill_manual(values = group_colors) +
  ylim(0, 8000) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "nFeature_RNA por grupo", y = "nFeature_RNA", x = "Grupo")


# --- nCOUNTS ---
vln_nCount <- ggplot(df_nCount, aes(x = amostra, y = nCount_RNA, fill = grupo)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_jitter(width = 0.2, size = 0.2, alpha = 0.1) +
  scale_fill_manual(values = group_colors) +
  ylim(0, 35000) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "nCount_RNA", y = " ", x = " ")

vln_nCount_group <- ggplot(df_nCount, aes(x = grupo, y = nCount_RNA, fill = grupo)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_jitter(width = 0.2, size = 0.2, alpha = 0.1) +
  scale_fill_manual(values = group_colors) +
  ylim(0, 35000) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "nCount_RNA por grupo", y = "nCount_RNA", x = "Grupo")


# ---- percent MT ---

vln_percentMT <- ggplot(df_percent, aes(x = amostra, y = percent.mt, fill = grupo)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_jitter(width = 0.2, size = 0.2, alpha = 0.1) +
  scale_fill_manual(values = group_colors) +
  ylim(0, 50) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Percent MT", y = NULL, x = NULL)

vln_percentMT_group <- ggplot(df_percent, aes(x = grupo, y = percent.mt, fill = grupo)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_jitter(width = 0.2, size = 0.2, alpha = 0.1) +
  scale_fill_manual(values = group_colors) +
  ylim(0, 50) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Percentual de genes mitocondriais por grupo", y = "percent.mt", x = "Grupo")

# ---- scater plot nCount vs nFeature ---

# --- todos integrado ---
scatter_qc <- ggplot(df_qc, aes(x = nCount_RNA, y = nFeature_RNA, color = grupo)) +
  geom_point(alpha = 0.3, size = 0.5) +
  scale_color_manual(values = group_colors) +
  facet_wrap(~amostra, scales = "free") +
  theme_minimal() +
  labs(title = "nFeature vs nCount por amostra", x = "nCount_RNA", y = "nFeature_RNA")

# ordena pelas cores dos grupos
ordem_amostras <- names(group_colors)

# cria coluna com os grupos
df_grupo <- df_qc[grepl(grupo_base, df_qc$grupo), ]

#lista de itens para o loop
grupos_desejados <- c("dt1", "dt2", "dt3", "dt4")

# --- Loop para gerar o gráfico e salvar ---
for (grupo_base in grupos_desejados) {
  
  df_grupo <- df_qc[grepl(grupo_base, df_qc$grupo), ]
  
  # Subgrupos daquele grupo (na ordem definida)
  subgrupos <- ordem_amostras[grepl(grupo_base, ordem_amostras)]
  
  # Reordenar os níveis da amostra com base nos subgrupos
  df_grupo$amostra <- factor(df_grupo$amostra, levels = unique(df_grupo$amostra[subgrupos %in% df_grupo$grupo]))
  
  # Gerar o gráfico
  p <- ggplot(df_grupo, aes(x = nCount_RNA, y = nFeature_RNA, color = grupo)) +
    geom_point(alpha = 0.3, size = 0.5) +
    scale_color_manual(values = group_colors) +
    facet_wrap(~amostra, scales = "free") +
    theme_minimal() +
    labs(title = paste("nFeature vs nCount (", grupo_base, ")", sep = ""),
         x = "nCount_RNA",
         y = "nFeature_RNA")
  
  print(p)
  
  ggsave(
    filename = file.path(output_dir, paste0("scatter_", grupo_base, "_por_amostra.png")),
    plot = p,
    width = 8,
    height = 6
  )
}


# -------------------------------
# 5. Salvar os plots
# -------------------------------

ggsave(file.path(output_dir, "all_violin_nFeature.png"), vln_nFeature, width = 10, height = 6, dpi = 300)
ggsave(file.path(output_dir, "all_violin_nCount.png"), vln_nCount, width = 10, height = 6, dpi = 300)
ggsave(file.path(output_dir, "all_violin_percentMT.png"), vln_percentMT, width = 10, height = 6, dpi = 300)
ggsave(file.path(output_dir, "all_scatter_qc.png"), scatter_qc, width = 14, height = 6, dpi = 300)

ggsave(file.path(output_dir, "group_violin_nFeature.png"), vln_nFeature_group, width = 10, height = 6, dpi = 300)
ggsave(file.path(output_dir, "group_violin_nCount.png"), vln_nCount_group, width = 10, height = 6, dpi = 300)
ggsave(file.path(output_dir, "group_violin_percentMT.png"), vln_percentMT_group, width = 10, height = 6, dpi = 300)


