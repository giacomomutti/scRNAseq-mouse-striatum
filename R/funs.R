colors <- c(
  "#FFC312", "#C4E538", "#12CBC4", "#FDA7DF", "#ED4C67",
  "#e17055", "#3c6382", "#e55039", "#0a3d62", "#d63031",
  "#EE5A24", "#009432", "#0652DD", "#9980FA", "#833471",
  "#EA2027", "#006266", "#1B1464", "#5758BB", "#6F1E51",
  "#40407a", "#706fd3", "#f7f1e3", "#34ace0", "#33d9b2",
  "#F79F1F", "#A3CB38", "#1289A7", "#D980FA", "#B53471",
  "#2c2c54", "#474787", "#aaa69d", "#227093", "#218c74",
  "#ff5252", "#ff793f", "#d1ccc0", "#ffb142", "#ffda79",
  "#b33939", "#cd6133", "#84817a", "#cc8e35", "#ccae62",
  "#00b894", "#ff7675", "#6c5ce7", "#636e72", "#fdcb6e"
)

colors_2 <- colors[c(1, 2, 3, 4, 5, 12, 17, 18, 21, 36, 50)]

vln_plot <- function(df, x, y,
                     threshold = quantile(df[, y], 0.1),
                     title = y,
                     quantiles = c(.1, .5, .9),
                     # group = NULL,
                     color = "#1F9A8AFF",
                     jitter = FALSE,
                     log = FALSE) {
  p <- df %>% ggplot(aes_string(x = x, y = y, fill = x, color = x)) +
    # this can be fixed, one condition more consequences. but now i don't have time.
    # computing quantiles gives a warning but no problem
    geom_violin(color = "black", draw_quantiles = quantiles, trim = FALSE, scale = "area", na.rm = TRUE) +
    geom_hline(yintercept = threshold, color = "firebrick") +
    {
      if (jitter) {
        geom_jitter(
          alpha = 0.3,
          size = 0.15
        )
      }
    } +
    {
      if (log) scale_y_log10(labels = scales::comma)
    } +
    # {if (log) annotation_logticks(sides = "l")} +
    {
      if (log == FALSE) scale_y_continuous()
    } +
    scale_x_discrete() +
    labs(title = title, x = "", y = "") +
    theme_bw() +
    theme(
      legend.position = "none"
    ) +
    coord_flip()
  p
}

# breaks = scales::trans_breaks("log10", function(x) 10^x, n = 4),
# labels = scales::trans_format("log10", scales::math_format(10^.x))

scatter_plot <- function(df, x, y, col,
                         v_threshold = quantile(df[, x], 0.1),
                         h_threshold = quantile(df[, y], 0.1),
                         legend_x = 0.9,
                         ...) {
  p <- df %>% ggplot(mapping = aes_string(x = x, y = y, color = col)) +
    geom_point(size = 0.2) +
    geom_vline(xintercept = v_threshold, color = "firebrick") +
    geom_hline(yintercept = h_threshold, color = "firebrick") +
    theme_bw() +
    ggpubr::stat_cor(aes(label = ..r.label..), method = "pearson", label.x.npc = "center") +
    theme(
      legend.position = c(legend_x, 0.5),
      legend.background = element_rect(fill = NA)
    )
}


par_coor <- function(df, subset = everything()) {
  df <- df %>%
    select(all_of(subset)) %>%
    select_if(is.numeric) %>%
    transmute_if(is.numeric, scales::rescale, to = c(0, 100))
  df %>%
    mutate(group = case_when(
      df[, 1] > quantile(df[, 1], 0.75) ~ "1 Qu.",
      df[, 1] > quantile(df[, 1], 0.50) &
        df[, 1] <= quantile(df[, 1], 0.75) ~ "2 Qu.",
      df[, 1] > quantile(df[, 1], 0.25) &
        df[, 1] <= quantile(df[, 1], 0.50) ~ "3 Qu.",
      df[, 1] <= quantile(df[, 1], 0.25) ~ "4 Qu."
    )) %>%
    mutate(ID = 1:n()) %>%
    pivot_longer(-c(ID, group), names_to = "key", values_to = "value") %>%
    ggplot(aes(y = value, x = key, group = ID, color = group)) +
    geom_line(alpha = 0.2, size = 0.1) +
    viridis::scale_color_viridis(
      discrete = TRUE,
      guide = guide_legend(override.aes = list(size = 1, alpha = 1))
    ) +
    scale_y_continuous(expand = c(0, 5)) +
    scale_x_discrete(expand = c(0.03, 0.03), limits = colnames(df)) +
    theme_bw() +
    labs(x = "", y = "Scaled Value") +
    theme(
      legend.title = element_blank(),
      legend.position = c(0.5, 0.98),
      legend.direction = "horizontal",
      legend.key = element_blank(),
      legend.background = element_rect(fill = NA)
    )
}

density_plot <- function(df, x, threshold = quantile(df[, x], 0.1),
                         color = "#1F9A8AFF",
                         title = x,
                         size = 1.5) {
  df %>% ggplot(aes_string(x = x)) +
    ggtitle(title) +
    geom_density(color = color, size = size) +
    theme_bw() +
    geom_vline(xintercept = threshold, size = size, color = "firebrick")
}

order_plot <- function(df, y, title = y,
                       threshold = quantile(df[, y], 0.1),
                       color = "#1F9A8AFF") {
  df %>%
    arrange(desc(df[, y])) %>%
    ggplot(aes_string(x = seq(1, length(df[, y])), y = y)) +
    geom_point(color = color) +
    # geom_point(aes_string( y = df2[,4]), color = 'purple', alpha = 0.1) +
    # geom_point(aes_string( y = df2[,3]), color = 'orange', alpha = 0.1) +
    geom_hline(yintercept = threshold, color = "firebrick") +
    labs(title = title, x = "Rank", y = "") +
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x, n = 4),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    theme_bw() +
    annotation_logticks(sides = "l")
}


prelim_analysis <- function(seurat_obj,
                            depth_threshold = NULL,
                            genes_threshold = NULL,
                            mit_threshold = NULL,
                            novelty_threshold = NULL,
                            main_color = "#1F9A8AFF") {
  df <- seurat_obj[[]]

  ifelse(is.null(genes_threshold), genes_threshold <- quantile(df$nFeature_RNA, 0.1), genes_threshold)
  ifelse(is.null(depth_threshold), depth_threshold <- quantile(df$nCount_RNA, 0.1), depth_threshold)
  ifelse(is.null(mit_threshold), mit_threshold <- quantile(df$percent.mt, 0.1), mit_threshold)
  ifelse(is.null(novelty_threshold), novelty_threshold <- quantile(df$novelty, 0.1), novelty_threshold)

  vln_theme <- list(
    scale_fill_manual(values = c(main_color)),
    scale_color_manual(values = c(main_color)),
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
  )

  p1.1 <- vln_plot(df, "orig.ident", "nCount_RNA", threshold = depth_threshold, log = TRUE, jitter = FALSE) +
    vln_theme
  p1.2 <- vln_plot(df, "orig.ident", "nFeature_RNA", threshold = genes_threshold, log = TRUE, jitter = FALSE) +
    vln_theme
  p1.3 <- vln_plot(df, "orig.ident", "percent.mt", threshold = mit_threshold, jitter = FALSE) +
    vln_theme

  table <- data.frame(signif(rbind(
    summary(df$nCount_RNA),
    summary(df$nFeature_RNA),
    summary(df$percent.mt),
    summary(df$novelty)
  ), 2))

  rownames(table) <- colnames(df)[2:5]
  colnames(table)[c(2, 5)] <- c("1st. Qu.", "2nd. Qu.")
  plot_tab <- ggpubr::ggtexttable(table)

  p2 <- par_coor(head(df, 1000))

  cont <- list(
    scale_color_viridis(discrete = FALSE),
    guides(color = guide_colourbar(
      barwidth = 0.2, barheight = 10,
      title.position = "left", title.theme = element_text(angle = 90, hjust = 0.5)
    ))
  )

  p3 <- scatter_plot(df, "nCount_RNA", "nFeature_RNA", "percent.mt",
    v_threshold = depth_threshold, h_threshold = genes_threshold
  ) + cont

  p3 <- ggMarginal(p3, type = "densigram", fill = main_color, alpha = 0.3)
  # p3.2 <- scatter_plot(df, 2,4,3, v_threshold = depth_threshold, h_threshold = mit_threshold)
  # p3.3 <- scatter_plot(df, 3,4,2, v_threshold = genes_threshold, h_threshold = mit_threshold)
  p4 <- order_plot(df, "nCount_RNA", color = main_color)

  p5 <- density_plot(df, "novelty", threshold = novelty_threshold, color = main_color)


  fst_row <- plot_grid(plot_grid(p1.1, p1.2, p1.3, ncol = 1), plot_tab, p2,
    ncol = 3,
    # widths = c(1, 1, 1, 3, 3),
    labels = c("A", "B", "C")
  )

  snd_row <- plot_grid(p3, p4, p5, nrow = 1, labels = c("D", "E", "F"))

  plot_grid(fst_row, snd_row, nrow = 2)
}

varfeat <- function(sm_object) {
  top10_strict <- head(VariableFeatures(sm_object), 10)

  # plot variable features with and without labels
  varfeat <- VariableFeaturePlot(sm_object, cols = c("#273c75", "#c23616"), pt.size = 0.5) +
    theme_bw() +
    labs(subtitle = "Most variable genes found with vst method") +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      plot.subtitle = element_text(face = "bold")
    )
  varfeat <- LabelPoints(plot = varfeat, points = top10_strict, repel = TRUE, xnudge = 0, ynudge = 0)
  varfeat
}


doublet_analysis <- function(seurat_df, col = c("firebrick", "navyblue")) {
  p1 <- seurat_df %>%
    group_by(multiplet_class) %>%
    summarise(count = n()) %>%
    mutate(percent = count / nrow(seurat_df)) %>%
    ggplot(aes(x = multiplet_class, fill = multiplet_class, y = count)) +
    geom_bar(stat = "identity", col = "black") +
    geom_text(aes(label = scales::percent(percent)), hjust = -.1, size = 4) +
    scale_fill_manual(values = col) +
    coord_flip() +
    labs(x = "", y = "", title = "Singlets vs doublets count") +
    theme_bw() +
    theme(
      legend.position = "none",
      title = element_text(face = "bold")
    )

  p2 <- vln_plot(seurat_df, "multiplet_class", "nCount_RNA", log = TRUE, title = "Number of transcripts") +
    scale_fill_manual(values = col) + labs(subtitle = "Log scale") +
    theme(plot.title = element_text(hjust = 1, face = "bold"))

  p3 <- vln_plot(seurat_df, "multiplet_class", "nCount_RNA", log = FALSE) +
    scale_fill_manual(values = col) + labs(subtitle = "Linear scale", title = "")

  p4 <- vln_plot(seurat_df, "multiplet_class", "nFeature_RNA", log = TRUE) +
    scale_fill_manual(values = col) + labs(subtitle = "Log scale", title = "Number of genes") +
    theme(plot.title = element_text(hjust = 1, face = "bold"))

  p5 <- vln_plot(seurat_df, "multiplet_class", "nFeature_RNA", log = FALSE) +
    theme(title = element_text(face = "bold")) +
    scale_fill_manual(values = col) + labs(subtitle = "Linear scale", title = "")

  p6 <- vln_plot(seurat_df, "multiplet_class", "percent.mt", log = FALSE, title = "Mito percentage") +
    scale_fill_manual(values = col) +
    theme(plot.title = element_text(face = "bold"))

  p7 <- vln_plot(seurat_df, "multiplet_class", "novelty", log = FALSE, title = "Novelty") +
    scale_fill_manual(values = col) +
    theme(plot.title = element_text(face = "bold"))

  bottom <- plot_grid(plot_grid(p2, p3, p4, p5, ncol = 4, nrow = 1),
    plot_grid(p6, p7),
    nrow = 2, ncol = 1, rel_heights = c(1, 1)
  )
  doublets_comparison <- plot_grid(p1, bottom, ncol = 1, rel_heights = c(0.3, 1))
  doublets_comparison
}


PC_elbow <- function(seurat, threshold = 10,
                     color = "firebrick", fill = "#273c75",
                     xmax = length(PCA_stdev)) {
  PCA_stdev <- seurat@reductions$pca@stdev
  xl <- as.character(seq(1, xmax))
  data.frame(
    PC = 1:length(PCA_stdev),
    stdev = PCA_stdev
  ) %>%
    ggplot(aes(x = factor(PC), y = stdev)) +
    geom_bar(
      stat = "identity", fill = fill,
      color = "black", alpha = 0.6
    ) +
    geom_vline(xintercept = threshold, color = color) +
    theme_bw() +
    xlim(xl) +
    labs(
      x = "Principal components",
      y = "Standard deviation",
      subtitle = "PC importance plot"
    ) +
    theme(
      plot.subtitle = element_text(face = "bold"),
      panel.grid.major.x = element_blank(),
      axis.text.x = element_text(angle = 90)
    )
}

# clustering functions

cluster_abundance <- function(seurat_object, res = FALSE) {
  df <- seurat_object[[]]
  if (res == TRUE) {
    df <- df %>%
      rownames_to_column(var = "cell") %>%
      pivot_longer(contains("snn"), names_to = "res", values_to = "cluster")
    df$res <- gsub("RNA_snn_res.", "", df$res)

    df %>%
      group_by(res, cluster) %>%
      tally() %>%
      ggplot(aes(x = cluster, y = n, fill = res)) +
      geom_bar(stat = "identity") +
      facet_grid(~res, scales = "free_x") +
      theme_bw() +
      theme(legend.position = "none")
  } else {
    df <- df
    df %>%
      group_by(seurat_clusters) %>%
      tally() %>%
      ggplot(aes(x = seurat_clusters, y = n, fill = seurat_clusters)) +
      geom_bar(stat = "identity") +
      theme_bw() +
      theme(legend.position = "none")
  }
}


cluster_density <- function(seurat_obj, y_id, res = FALSE) {
  df <- seurat_obj[[]]

  if (res == TRUE) {
    
    df <- df %>%
      rownames_to_column(var = "cell") %>%
      pivot_longer(contains("snn"), names_to = "res", values_to = "cluster")
    df$res <- gsub("RNA_snn_res.", "", df$res)

    ggplot() +
      geom_jitter(
        data = df, aes(x = cluster, y = .data[[y_id]], color = cluster, fill = cluster),
        alpha = 0.25, size = 0.1, show.legend = FALSE
      ) +
      geom_half_violin(
        data = df, aes(cluster, y = .data[[y_id]], fill = cluster),
        side = "l", alpha = 0.4, show.legend = FALSE, scale = "width"
      ) +
      geom_half_boxplot(
        data = df, aes(cluster, y = .data[[y_id]], fill = cluster),
        side = "r", alpha = 0.4, show.legend = FALSE
      ) +
      facet_grid(rows = vars(res)) +
      theme_bw() +
      labs(x = "", y = "Number of expressed genes") +
      theme(panel.grid.major.x = element_blank())
    
  } else {
    df <- df
    df %>% ggplot() +
      geom_jitter(
        data = df, aes(
          x = seurat_clusters, y = .data[[y_id]],
          color = seurat_clusters, fill = seurat_clusters
        ),
        alpha = 0.25, size = 0.1, show.legend = FALSE
      ) +
      geom_half_violin(
        data = df, aes(seurat_clusters, y = .data[[y_id]], fill = seurat_clusters),
        side = "l", alpha = 0.4, show.legend = FALSE, scale = "width"
      ) +
      geom_half_boxplot(
        data = df, aes(seurat_clusters, y = .data[[y_id]], fill = seurat_clusters),
        side = "r", alpha = 0.4, show.legend = FALSE
      ) +
      theme_bw() +
      labs(x = "", y = "Number of expressed genes") +
      theme(panel.grid.major.x = element_blank())
  }
}

# TODO venndiagram like graph
silhouette_analysis <- function(sm_obj, npcs = 15) {
  df <- sm_obj@meta.data
  # silhouette
  # https://github.com/satijalab/Integration2019/blob/master/analysis_code/integration/integration_metrics.R#L36

  distance_matrix <- dist(Embeddings(sm_obj[["pca"]])[, 1:npcs])

  for (col in grep("snn", colnames(df))) {
    res <- paste(sub(".*?\\.", "", colnames(df)[col]))
    clusters <- df[, col]
    silhouette <- silhouette(as.numeric(clusters), dist = distance_matrix)
    df[, paste0("ss", res)] <- silhouette[, 3]
  }


  plot_list <- list()
  colnm <- colnames(df)
  num <- 1

  for (i in grep("^ss", colnm)) {
    col <- colnm[i]
    res <- gsub("[^0-9.]", "", col)
    cluster <- colnm[grep(res, colnm)][1]

    mean_silhouette_score <- mean(df[, i])
    p <- df %>%
      mutate(barcode = rownames(.)) %>%
      arrange_(cluster, col) %>%
      mutate(barcode = factor(barcode, levels = barcode)) %>%
      ggplot() +
      geom_col(aes_string("barcode", col, fill = cluster), show.legend = FALSE) +
      geom_hline(yintercept = mean_silhouette_score, color = "red", linetype = "dashed") +
      scale_x_discrete(name = "Cells") +
      labs(title = paste("Resolution: ", res)) +
      scale_y_continuous(name = "Silhouette score") +
      scale_fill_manual(values = colors) +
      theme_bw() +
      theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
    plot_list[[num]] <- p
    num <- num + 1
  }

  plot_grid(plotlist = plot_list, nrow = 1)
}


# this function will be really unstable, column must have same name
# don't have time to make it stable

dimplot_custom <- function(seurat_obj, color, red, idx = c(1, 2)) {
  coor <- seurat_obj[[red]]@cell.embeddings[, idx]
  colnm <- colnames(coor)
  df <- cbind(seurat_obj@meta.data, coor)
  df %>% ggplot(aes_string(x = colnm[1], y = colnm[2], color = color)) +
    geom_point(size = 0.2) +
    theme_bw() +
    coord_fixed() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      aspect.ratio = 1
    )
}


celltype_dimred <- function(seurat_obj, group, red, idx = c(1, 2), facet = FALSE, annotate = TRUE) {
  centers <- tibble(
    fst = seurat_obj[[red]]@cell.embeddings[, idx[1]],
    snd = seurat_obj[[red]]@cell.embeddings[, idx[2]],
    cell_type = seurat_obj[[]][, group]
  ) %>%
    group_by(cell_type) %>%
    summarize(x = median(fst), y = median(snd))

  dimplot_custom(seurat_obj, group, red = red, idx = idx) +
    {
      if (facet) facet_grid(as.formula(paste("~", group)))
    } +
    {
      if (annotate) {
        ggrepel::geom_label_repel(
          data = centers,
          mapping = aes(x, y, label = cell_type),
          size = 2,
          fill = "white",
          color = "black",
          fontface = "bold",
          alpha = 0.7,
          show.legend = FALSE
        )
      }
    } +
    guides(col = guide_legend(override.aes = list(size = 4)))
}


color_phase <- colors[c(5, 26, 28, 32)]

phase_prop <- function(seurat_obj, by.cluster = FALSE) {
  df <- seurat_obj@meta.data

  if (by.cluster) {
    cluster_bycycle <- df %>%
      group_by(seurat_clusters, Phase) %>%
      tally() %>%
      ungroup()

    cluster_bycycle %>% ggplot(aes(seurat_clusters, n, fill = Phase)) +
      geom_bar(stat = "identity", position = "stack") +
      scale_fill_manual(values = color_phase) +
      scale_y_continuous(expand = c(0.01, 0.01)) +
      theme_bw() +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
  } else {
    df %>%
      group_by(Phase) %>%
      tally() %>%
      mutate(proj = "sc") %>%
      ggplot(aes(x = proj, y = n, fill = Phase)) +
      geom_bar(position = "stack", stat = "identity") +
      geom_text(aes(label = Phase), position = position_stack(vjust = 0.5), size = 3, color = "grey18") +
      scale_fill_manual(values = color_phase) +
      theme_bw() +
      theme(
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
  }
}


dimred_analysis <- function(sm_object) {
  plot_list <- list()

  cont <- list(
    scale_color_viridis(
      guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"),
      labels = scales::comma
    ),
    theme(legend.position = "left")
  )

  disc <- list(
    guides(col = guide_legend(override.aes = list(size = 2))),
    theme(legend.position = "left")
  )

  u1 <- dimplot_custom(sm_object, "nFeature_RNA", "umap") +
    labs(color = "Number of\ntranscripts:") +
    cont

  t1 <- dimplot_custom(sm_object, "nFeature_RNA", "tsne") +
    cont +
    theme(legend.position = "none")

  pc1 <- dimplot_custom(sm_object, "nFeature_RNA", "pca") +
    cont +
    theme(legend.position = "none")

  row1 <- u1 | t1 | pc1

  u2 <- dimplot_custom(sm_object, "percent.mt", "umap") +
    labs(color = "Mito %:") +
    cont

  t2 <- dimplot_custom(sm_object, "percent.mt", "tsne") +
    cont +
    theme(legend.position = "none")

  pc2 <- dimplot_custom(sm_object, "percent.mt", "pca") +
    cont +
    theme(legend.position = "none")

  row2 <- u2 | t2 | pc2

  u3 <- dimplot_custom(sm_object, "novelty", "umap") +
    labs(color = "Novelty:") +
    cont

  t3 <- dimplot_custom(sm_object, "novelty", "tsne") +
    cont +
    theme(legend.position = "none")

  pc3 <- dimplot_custom(sm_object, "novelty", "pca") +
    cont +
    theme(legend.position = "none")

  row3 <- u3 | t3 | pc3

  plot_list[["cont"]] <- plot_grid(row1, row2, row3, nrow = 3, align = "v")

  phase_u <- dimplot_custom(sm_object, "Phase", "umap") +
    labs(color = "Phase:") +
    scale_color_manual(values = color_phase) +
    disc
  phase_t <- dimplot_custom(sm_object, "Phase", "tsne") +
    labs(color = "Phase:") +
    scale_color_manual(values = color_phase) +
    disc +
    theme(legend.position = "none")
  phase_pc <- dimplot_custom(sm_object, "Phase", "pca") +
    labs(color = "Phase:") +
    scale_color_manual(values = color_phase) +
    disc +
    theme(legend.position = "none")

  phase <- phase_u | phase_t | phase_pc

  plot_list[["phase"]] <- phase

  phase2 <- (phase_t + facet_grid(~Phase)+theme(legend.position='none')) / 
    (phase_u + facet_grid(~Phase)+theme(legend.position='none')) / 
    (phase_pc + facet_grid(~Phase)+theme(legend.position='none'))

  plot_list[["phase_facet"]] <- phase2


  clu_u <- dimplot_custom(sm_object, "seurat_clusters", "umap") +
    labs(color = "Cluster: ") +
    scale_color_manual(values = colors) +
    disc

  clu_t <- dimplot_custom(sm_object, "seurat_clusters", "tsne") +
    labs(color = "Cluster: ") +
    scale_color_manual(values = colors) +
    disc +
    theme(legend.position = "none")

  clu_pc <- dimplot_custom(sm_object, "seurat_clusters", "pca") +
    labs(color = "Cluster: ") +
    scale_color_manual(values = colors) +
    disc +
    theme(legend.position = "none")

  clu <- clu_u | clu_t | clu_pc

  plot_list[["clust"]] <- clu

  # res_list <- list()
  # 
  # colnm <- colnames(sm_object@meta.data)
  # num <- 1
  # 
  # for (i in grep("snn_res", colnm)) {
  #   col <- colnm[i]
  #   res <- gsub("[^0-9.]", "", col)
  #   res <- sub(".", "", res)
  # 
  #   pu <- dimplot_custom(sm_object, col, "umap") +
  #     labs(color = "Cluster: ") +
  #     scale_color_manual(values = colors) +
  #     disc +
  #     ggtitle(paste("Resolution: ", res))
  #   pt <- dimplot_custom(sm_object, col, "tsne") +
  #     labs(color = "Cluster: ") +
  #     scale_color_manual(values = colors) +
  #     disc +
  #     ggtitle(paste("Resolution: ", res)) +
  #     theme(legend.position = "none")
  #   pc <- dimplot_custom(sm_object, col, "pca") +
  #     labs(color = "Cluster: ") +
  #     scale_color_manual(values = colors) +
  #     disc +
  #     ggtitle(paste("Resolution: ", res)) +
  #     theme(legend.position = "none")
  # 
  #   row <- pu | pt | pc
  # 
  #   res_list[[num]] <- row
  #   num <- num + 1
  # }
  # 
  # res <- plot_grid(plotlist = res_list, nrow = num - 1, align = "v")
  # 
  # plot_list[["res"]] <- res
  # 
  invisible(plot_list)
}

top_spread <- function(markers, n) {
  markers %>%
    group_by(cluster) %>%
    top_n(n = n, wt = avg_logFC) %>%
    mutate(rank = rank(-avg_logFC, ties.method = "first")) %>%
    ungroup() %>%
    select(rank, cluster, gene) %>%
    spread(key = rank, value = gene) %>%
    rename_at(vars(-1), ~ paste0("rank", .))
}

panglao_searcher <- function(top_output, collapse = " ", append = "type") {
  df <- top_output
  clust <- 1
  colidx <- grep("rank", colnames(df))
  df[append] <- NA
  df <- df %>% mutate_if(is.logical, as.character)
  # useful shortcuts


  while (clust <= nrow(df)) {
    # ol <- 'Oligodendrocytes'
    # n <- 'Neurons'
    # u <- 'Unknown'
    # a <- 'Astrocytes'
    # m <- 'Microglia'
    # e <- 'Endothelial cells'
    # opc <- 'Oligodendrocytes parental cells'
    print(paste0("Searching cell type of cluster:", clust))

    mark <- paste0(as.character(unlist(df[clust, colidx])), collapse = collapse)

    # write to clipboard to search in panglao
    writeClipboard(mark)

    print("Top 3 markers copied to clipboard, search on panglao and input class:")
    var <- readline()

    var <- get(var)

    df[clust, append] <- var
    clust <- +clust + 1
  }
  return(df)
}


celltype_analysis <- function(seurat_obj, group) {
  plot_list <- list()

  co <- list(
    scale_color_manual(values = colors),
    theme(legend.position = "none")
  )
  
  p1 <- celltype_dimred(seurat_obj, group, red = "umap") + co
  p2 <- celltype_dimred(seurat_obj, group, red = "tsne") + co
  p3 <- celltype_dimred(seurat_obj, group, red = "pca") + co 
  
  fst <- p1 | p2 | p3
  
  plot_list[["clust"]] <- fst

  p4 <- celltype_dimred(seurat_obj, group, red = "umap", annotate = F, facet = T) + co
  p5 <- celltype_dimred(seurat_obj, group, red = "tsne", annotate = F, facet = T) + co
  p6 <- celltype_dimred(seurat_obj, group, red = "pca", annotate = F, facet = T) + co

  snd <- p4 / p5 / p6
  
  plot_list[["clust_facet"]] <- snd
  
  invisible(plot_list)
}
