# 分步工作流 -----------

# 从当前 R 环境中移除（卸载）已经加载的包
detach("package:ARBO", unload = TRUE, character.only = TRUE)

# 强制重新加载
rm(list = ls())
library(devtools)
# 加载你正在开发的本地 R 包源码，并模拟包安装后的状态
devtools::load_all("/home/hsinyinteng/ARBO", reset = TRUE)


library(Cardinal)
Eye <- readImzML("/home/hsinyinteng/Spatial_Metabolomics/fish_eye/01_root_mean_square.imzML")
Eye <- Eye |> peakPick(SNR = 2) |> peakAlign()
Eye
image(Eye)

# 1. 提取矩阵
dat <- extract_spectra_matrix(Eye)
print(dim(dat$spectra))
print(head(dat$pixel_info))

# 2. 数据检查
check_spectra_matrix(dat$spectra)

# 3. 去常数列
spectra_filtered <- remove_constant_features(dat$spectra)
print(dim(spectra_filtered))

# 4. Min-max归一化
spectra_scaled <- apply_feature_scaling(spectra_filtered, method = "minmax")
print(range(spectra_scaled, na.rm = TRUE))

# 5. 检查Python环境
check_umap_python_env("/home/hsinyinteng/miniconda3/envs/dxy_python9/bin/python")

# 6. UMAP
umap_df <- run_umap_py(
  x = spectra_scaled,
  python_path = "/home/hsinyinteng/miniconda3/envs/dxy_python9/bin/python",
  metric = "cosine",
  n_neighbors = 10L,
  min_dist = 0.05,
  n_components = 2L,
  random_state = 2025L
)

print(head(umap_df))
print(dim(umap_df))

# 7. kmeans
km_res <- run_clustering(
  embedding = umap_df,
  method = "kmeans",
  centers = 10
)
print(table(km_res$cluster))

# 8. 合并结果
cluster_df <- build_cluster_dataframe(
  embedding = umap_df,
  cluster = km_res$cluster,
  pixel_info = dat$pixel_info
)

print(head(cluster_df))
cluster_df$kmeans_cluster <- cluster_df$cluster

ggplot(cluster_df, aes(x = x, y = y, color = factor(kmeans_cluster))) +
  geom_point(size = 1, alpha = 0.9) +
  scale_y_reverse() +
  scale_color_manual(
    values = c(
      "1" = "#8B0000", "2" = "#2E8B57", "3" = "#FF0000",
      "4" = "#FFC0CB", "5" = "#FF00FF", "6" = "#00FFFF",
      "7" = "#A9A9A9", "8" = "#FF7F0E", "9" = "#9370DB",
      "10" = "#FFD700", "11" = "#00FF00", "12" = "#000000",
      "13" = "#0000FF"
    ),
    name = "Kmeans Cluster",
    guide = guide_legend(
      title.position = "left",       # 标题位置
      title.hjust = 0.5,             # 标题居中
      byrow = TRUE,                  # 按行填充
      override.aes = list(size = 8), # 增大图例点大小
      ncol=4                         # 图例的列数
    )
  ) +
  coord_fixed(ratio = 1) +
  theme_minimal(base_family = "Arial") +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "bottom",    # 图例放在底部
    legend.box = "horizontal",     # 水平排列
    legend.justification = "center",
    legend.spacing.x = unit(0.5, 'cm'), # 图例子项水平间距
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 20),
    plot.margin = margin(10, 10, 20, 10) # 调整底部边距避免图例被裁剪
  )


library(ggplot2)
ggplot(cluster_df, aes(x = UMAP1, y = UMAP2, color = kmeans_cluster)) +
  geom_point(size = 0.1, alpha = 0.7) +  # 使用点形状
  labs(# title = "UMAP Clustering of All Pixels",
    x = "UMAP1",
    y = "UMAP2") +
  scale_color_manual(
    values = c("1"  = "#8B0000","2"  = "#2E8B57","3"  = "#FF0000",
               "4"  = "#FFC0CB","5"  = "#FF00FF","6"  = "#00FFFF",
               "7"  = "#A9A9A9","8"  = "#FF7F0E","9"  = "#9370DB",
               "10" = "#FFD700","11" = "#00FF00","12" = "#000000",
               "13" = "#0000FF"),
    guide = guide_legend(override.aes = list(size = 8),# 调整图例点大小
                         ncol=1) # 图例的列数
  ) +
  theme_minimal(base_family = "Arial") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),  # 添加黑色边框
        panel.grid.major = element_blank(),  # 去除主要网格线
        panel.grid.minor = element_blank(),   # 去除次要网格线
        plot.title = element_text(hjust = 0.5, size = 16, face = "plain", family = "Arial"),
        legend.text = element_text(size = 18, face = "plain", family = "Arial"),
        legend.title = element_text(size = 20, face = "plain", family = "Arial"),
        axis.title.x = element_text(size = 18, face = "plain", family = "Arial"),
        axis.title.y = element_text(size = 18, face = "plain", family = "Arial"),
        axis.text = element_text(size = 16, face = "plain", family = "Arial")
  ) +
  coord_cartesian()  # 设置坐标轴范围




Eye2 <- attach_cluster_to_pixeldata(Eye, cluster_df)

head(as.data.frame(Cardinal::pixelData(Eye2)))
colnames(Cardinal::pixelData(Eye2))
head(as.data.frame(Cardinal::pixelData(Eye2))$pixel_ID)


# 一键工作流 --------------------
devtools::load_all("/home/hsinyinteng/ARBO", reset = TRUE)

library(Cardinal)
Eye <- readImzML("/home/hsinyinteng/Spatial_Metabolomics/fish_eye/01_root_mean_square.imzML")
Eye <- Eye |> peakPick(SNR = 2) |> peakAlign()

res <- spatial_clustering_workflow(
  msi_obj = Eye,
  python_path = "/home/hsinyinteng/miniconda3/envs/dxy_python9/bin/python",
  centers = 10
)

head(res$cluster_df)
head(as.data.frame(Cardinal::pixelData(res$msi_obj)))

cluster_df <- res$cluster_df
# umap_df <- res$umap_df
Eye <- res$msi_obj
cluster_df$kmeans_cluster <- cluster_df$cluster
pixelData(Eye)

ggplot(cluster_df, aes(x, y, color = factor(kmeans_cluster))) +
  geom_point(size = 1) +
  scale_y_reverse() +
  coord_fixed() +
  theme_void()

ggplot(cluster_df, aes(UMAP1, UMAP2, color = factor(kmeans_cluster))) +
  geom_point(size = 0.2) +
  theme_classic()

image(Eye)


# Spatially Enriched Metabolites(SEMs) screen ---------
devtools::load_all("/home/hsinyinteng/ARBO", reset = TRUE)

library(Cardinal)
data(cherry_tomato_msi)
tomato <- cherry_tomato_msi
tomato

res <- SEMs_screen(
  msi_obj = tomato,
  group_col = "cluster",
  t_statistic_threshold = 20,
  cor_threshold = 0.5,
  M1_threshold = 0.5,
  M2_threshold = 0.5,
  save_table = TRUE,
  file_name = "Spatially_Enriched_Metabolites",
  file_format = "xlsx"
)

dim(res$results_df)
dim(res$SEMs_df)

image(tomato,
      mz=sort(res$SEMs_df$mz),
      smooth="gaussian",
      enhance="hist",
      scale=T,
      key=F)

# SEMs 重构 -----------
tomato_SEMs_only <- subsetFeatures(tomato,
                                   mz %in% sort(res$SEMs_df$mz))

res <- spatial_clustering_workflow(
  msi_obj = tomato_SEMs_only,
  python_path = "/home/hsinyinteng/miniconda3/envs/dxy_python9/bin/python",
  centers = 2
)
cluster_df <- res$cluster_df

library(ggplot2)
ggplot(cluster_df, aes(x, y, color = factor(cluster))) +
  geom_point(size = 0.1) +
  scale_y_reverse() +
  coord_fixed() +
  theme_void()

# image overlay -------
devtools::load_all("/home/hsinyinteng/ARBO", reset = TRUE)

library(Cardinal)
data(cherry_tomato_msi)
tomato <- cherry_tomato_msi
tomato
image(tomato, mz=featureData(tomato)$mz)
tomato <- summarizePixels(tomato, c(TIC="sum"))

fruit <- subsetPixels(tomato, pixelData(tomato)$cluster == "Fruit")
fruit

msi_img_overlay(
  msi_data01 = tomato,
  value01 = "TIC",
  color01 = colorRampPalette(c("white", "grey"))(100),
  msi_data02 = fruit,
  value02 = 240.8123,
  color02 = colorRampPalette(c("black", "blue", "cyan", "yellow", "red"))(100),
  show_axes = FALSE
  )

msi_img_overlay(
  msi_data01 = tomato,
  value01 = "TIC",
  color01 = colorRampPalette(c("white", "grey"))(100),
  msi_data02 = fruit,
  value02 = "cluster",
  color02 = c(
    Calyx = "#E64B35",
    Fruit = "#3C5488"
  ),
  show_axes = FALSE
)

