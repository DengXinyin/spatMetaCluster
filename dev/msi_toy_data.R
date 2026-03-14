#' Toy spatial metabolomics dataset
#'
#' A simulated \code{MSImagingExperiment} object for demonstrating the
#' workflow in \pkg{spatMetaCluster}. The dataset contains a small 2D spatial
#' grid with region-specific spectral patterns.
#'
#' This dataset is intended for package examples, tutorials, and testing only.
#' It is not a real biological dataset and should not be used for biological
#' interpretation.
#'
#' @format An \code{MSImagingExperiment} object with spatial coordinates,
#' simulated m/z features, and a known region label stored in pixel metadata.
#' @source Simulated using \code{Cardinal::simulateSpectra()}.
"toy_msi"


# image process -------
library(dplyr)
library(ggplot2)
library(Cardinal)
library(jpeg)
# img <- readJPEG("/home/hsinyinteng/spatMetaCluster/tests/DSC00349.jpg")

library(png)
img <- readPNG("/home/hsinyinteng/spatMetaCluster/tests/Cherry_tomato.png")


# 像素降采样
step <- 4
img <- img[
  seq(1, dim(img)[1], by = step),
  seq(1, dim(img)[2], by = step),
  ,
  drop = FALSE
]


h <- dim(img)[1]
w <- dim(img)[2]
nch <- dim(img)[3]

pixel_df <- expand.grid(
  y = 1:h,
  x = 1:w
)

pixel_df$R01 <- img[cbind(pixel_df$y, pixel_df$x, 1)]
pixel_df$G01 <- img[cbind(pixel_df$y, pixel_df$x, 2)]
pixel_df$B01 <- img[cbind(pixel_df$y, pixel_df$x, 3)]

if (nch >= 4) {
  pixel_df$A01 <- img[cbind(pixel_df$y, pixel_df$x, 4)]
} else {
  pixel_df$A01 <- 1
}

pixel_df <- pixel_df %>%
  mutate(
    R = round(R01 * 255),
    G = round(G01 * 255),
    B = round(B01 * 255),
    A = A01
  ) %>%
  filter(
    A > 0.01,
    !(R >= 245 & G >= 245 & B >= 245)
  )


set.seed(123)
rgb_data <- pixel_df %>% select(R, G, B)
km <- kmeans(rgb_data, centers = 4, iter.max = 50)

pixel_df$cluster <- as.factor(km$cluster)

centers <- round(km$centers)
cluster_colors <- rgb(
  centers[,1], centers[,2], centers[,3],
  maxColorValue = 255
)

p1 <- ggplot(pixel_df, aes(x = x, y = y, fill = cluster)) +
  geom_raster() +
  scale_y_reverse() +
  coord_fixed() +
  theme_void() +
  scale_fill_manual(values = cluster_colors) +
  labs(fill = "Cluster")

print(p1)


# 生成 spectra ----------
generate_spectra_from_clusters <- function(
    pixel_df,
    n_features = 300,
    seed = 2025,
    base_mean = 200,
    base_sd = 30,
    cluster_boost = 800,
    noise_sd = 50,
    rgb_weight = 0.3,
    mz_min = 100,
    mz_max = 1000,
    min_cluster_gap = 100,
    peak_width = 8
) {
  set.seed(seed)

  cl_factor <- as.factor(pixel_df$cluster)
  cl <- as.integer(cl_factor)
  cluster_levels <- levels(cl_factor)
  k <- length(cluster_levels)
  n_pixels <- nrow(pixel_df)

  # 检查 m/z 范围是否足够容纳各 cluster 的最小间隔
  if ((k - 1) * min_cluster_gap > (mz_max - mz_min)) {
    stop("Cluster number is too large for the required min_cluster_gap within mz range.")
  }

  # m/z 轴
  mz <- seq(mz_min, mz_max, length.out = n_features)

  # 基础谱矩阵：pixel x feature
  spectra <- matrix(
    rnorm(n_pixels * n_features, mean = base_mean, sd = base_sd),
    nrow = n_pixels,
    ncol = n_features
  )
  spectra[spectra < 0] <- 0

  # 1. 给每个 cluster 分配一个特异峰中心，彼此至少差 min_cluster_gap
  cluster_centers_mz <- seq(
    from = mz_min + 50,
    by = min_cluster_gap,
    length.out = k
  )

  # 若最后超过 mz_max - 50，则在范围内均匀铺开
  if (max(cluster_centers_mz) > (mz_max - 50)) {
    cluster_centers_mz <- seq(
      from = mz_min + 50,
      to = mz_max - 50,
      length.out = k
    )
    # 再检查相邻间隔是否满足要求
    if (any(diff(cluster_centers_mz) < min_cluster_gap)) {
      stop("Cannot assign cluster-specific mz centers with the required min_cluster_gap.")
    }
  }

  names(cluster_centers_mz) <- cluster_levels

  # 找每个 cluster 对应中心附近的 feature 索引
  cluster_feature_blocks <- lapply(cluster_centers_mz, function(center_mz) {
    which(mz >= (center_mz - peak_width) & mz <= (center_mz + peak_width))
  })

  # 2. 加入 cluster 特异峰
  for (lv in cluster_levels) {
    idx_pixel <- which(cl_factor == lv)
    idx_feat <- cluster_feature_blocks[[lv]]

    if (length(idx_feat) > 0 && length(idx_pixel) > 0) {
      spectra[idx_pixel, idx_feat] <-
        spectra[idx_pixel, idx_feat] +
        matrix(
          rnorm(length(idx_pixel) * length(idx_feat),
                mean = cluster_boost,
                sd = cluster_boost * 0.15),
          nrow = length(idx_pixel),
          ncol = length(idx_feat)
        )
    }
  }

  # 3. RGB 调制
  R_scaled <- pixel_df$R / 255
  G_scaled <- pixel_df$G / 255
  B_scaled <- pixel_df$B / 255

  idx_r <- seq(1, n_features, by = 3)
  idx_g <- seq(2, n_features, by = 3)
  idx_b <- seq(3, n_features, by = 3)

  spectra[, idx_r] <- spectra[, idx_r] +
    outer(R_scaled, rep(rgb_weight * 300, length(idx_r)))
  spectra[, idx_g] <- spectra[, idx_g] +
    outer(G_scaled, rep(rgb_weight * 300, length(idx_g)))
  spectra[, idx_b] <- spectra[, idx_b] +
    outer(B_scaled, rep(rgb_weight * 300, length(idx_b)))

  # 4. 加入噪声
  spectra <- spectra +
    matrix(
      rnorm(n_pixels * n_features, mean = 0, sd = noise_sd),
      nrow = n_pixels,
      ncol = n_features
    )

  spectra[spectra < 0] <- 0

  colnames(spectra) <- paste0("mz_", sprintf("%.2f", mz))
  rownames(spectra) <- paste0("pixel_", seq_len(n_pixels))

  list(
    spectra = spectra,
    mz = mz,
    cluster_centers_mz = cluster_centers_mz,
    cluster_feature_blocks = cluster_feature_blocks
  )
}

sim_res <- generate_spectra_from_clusters(
  pixel_df = pixel_df,
  n_features = 1000,
  seed = 2025,
  mz_min = 100,
  mz_max = 1000,
  min_cluster_gap = 100,
  peak_width = 10
)

spectra <- sim_res$spectra
mz <- sim_res$mz
sim_res$cluster_centers_mz



# 构建 MSI 对象 ------
make_demo_msi <- function(pixel_df, spectra, mz, run_name = "MSI_demo") {
  coord <- pixel_df[, c("x", "y")]
  run <- factor(rep(run_name, nrow(pixel_df)))

  fdata <- Cardinal::MassDataFrame(mz = mz)

  pdata <- Cardinal::PositionDataFrame(
    run = run,
    coord = coord
  )

  pdata$cluster <- pixel_df$cluster
  pdata$R <- pixel_df$R
  pdata$G <- pixel_df$G
  pdata$B <- pixel_df$B
  pdata$pixel_ID <- paste0(run_name, "_", seq_len(nrow(pixel_df)))

  Cardinal::MSImagingExperiment(
    spectraData = t(spectra),
    featureData = fdata,
    pixelData = pdata
  )
}

demo_msi <- make_demo_msi(pixel_df, spectra, mz, "cherry_tomato")
demo_msi

demo_msi2 <- demo_msi |> peakPick(SNR = 2) |> peakAlign()

demo_msi2 <- subsetFeatures(demo_msi2,
                            mz > 141 & mz < 160 |
                              mz > 240 & mz < 260 |
                              mz > 343 & mz < 358 |
                              mz > 443 & mz < 458)
demo_msi2

# cherry_tomato ---------
devtools::load_all("/home/hsinyinteng/spatMetaCluster", reset = TRUE)
data(cherry_tomato_msi)
demo_msi2 <- cherry_tomato_msi
demo_msi2


plot(demo_msi2, i=1020,
     isPeaks = T,
     annPeaks = 6)

image(demo_msi2,
      mz=featureData(demo_msi2)$mz,
      # mz=c(
      #   141.1193, 143.5521, 146.3088, 149.3382, 152.0502,
      #   154.7453, 157.0725, 159.1878, 241.2659, 244.0796,
      #   247.6973, 252.2827, 256.8508, 259.4595, 343.8279,
      #   351.3438, 357.4678, 444.2957, 451.9370, 457.5927),
      smooth="gaussian",
      enhance="hist",
      scale=T,
      key=F)


# 聚类方法对比 ---------
demo_msi2_pca <- PCA(demo_msi2, ncomp=4)
demo_msi2_pca
image(demo_msi2_pca, smooth="adaptive", enhance="histogram")


demo_msi2_nmf <- NMF(demo_msi2, ncomp=4, niter=30)
demo_msi2_nmf
image(demo_msi2_nmf, smooth="adaptive", enhance="histogram")


set.seed(2026)
demo_msi2_ssc <- spatialShrunkenCentroids(
  demo_msi2,
  weights="adaptive", r=2, k=4, s=2)

demo_msi2_ssc
image(demo_msi2_ssc, i=1:1)



devtools::load_all("/home/hsinyinteng/spatMetaCluster", reset = TRUE)
res <- spatial_clustering_workflow(
  msi_obj = demo_msi2,
  python_path = "/home/hsinyinteng/miniconda3/envs/dxy_python9/bin/python",
  centers = 4,
  metric = "cosine",
  n_neighbors = 25L,
  min_dist = 0.3,
  n_components = 2L
)
head(res$cluster_df)

cluster_df <- res$cluster_df
# umap_df <- res$umap_df
# duck_msi3 <- res$msi_obj
ggplot(cluster_df, aes(x, y, color = factor(kmeans_cluster))) +
  geom_point(size = 0.1) +
  scale_y_reverse() +
  coord_fixed() +
  theme_void() +
  labs(color = "Cluster") +
  guides(
    color = guide_legend(
      override.aes = list(size = 4)  # 图例点的大小
    )
  )

ggplot(cluster_df, aes(UMAP1, UMAP2, color = factor(kmeans_cluster))) +
  geom_point(size = 0.2) +
  theme_classic() +
  labs(color = "Cluster") +
  guides(
    color = guide_legend(
      override.aes = list(size = 4)  # 图例点的大小
    )
  )


ggplot(pixel_df, aes(x = x, y = y, fill = cluster)) +
  geom_raster() +
  scale_y_reverse() +
  coord_fixed() +
  theme_void() +
  scale_fill_manual(values = cluster_colors) +
  labs(fill = "Cluster")




cherry_tomato_msi <- demo_msi2
# 文件存放在 data/cherry_tomato_msi.rda，用户可使用data(demo_msi)进行调用
usethis::use_data(cherry_tomato_msi, overwrite = TRUE)
# 别忘了重新生成文档
# 如果你用了 roxygen2：
devtools::document()
# 检查帮助文件：
?cherry_tomato_msi
# 加载数据测试：
data(cherry_tomato_msi)
cherry_tomato_msi


# 保存结果至mouse_14_ovary_pixel_merged_1236_features.imzML文件夹，包括log文件, metadata文件, pdata文件，imzML文件和ibd文件
imzfile <- file.path("/home/hsinyinteng/Spatial_Metabolomics/mouse_ovary/continuous_result_20241113/", "mouse_14_ovary_pixel_merged_1236_features.imzML")
writeMSIData(demo_msi2, file = imzfile)

list.files(imzfile)
file.exists(imzfile)  # 返回 TRUE 表示文件生成成功
file.size(imzfile) > 0  # 返回 TRUE 表示文件非空

