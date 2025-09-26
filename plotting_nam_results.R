plotting_nam_results <- function(rils_dt = NULL, parents_dt = NULL, founders_dt = NULL,
                                 donor_parent_col = "Donor_parent", recur_parent_col = "Recurrent_parent",
                                 trait_name = "Trait1", trait_unit = "cm", n_plot_rows = 5, n_plot_cols = 2){
  ### You only need to provide one of parents_dt OR founders_dt. Not both.
  require(tidyverse)
  
  nam_dt <- rils_dt |> 
    dplyr::mutate(family = paste(!!sym(recur_parent_col), !!sym(donor_parent_col), sep = "X"),
                  gv = Trait1_GV, pheno = Trait1_Pheno, type = "RIL") |> 
    select(family, gv, pheno, type)
  
  nam_dt_summarized <- group_by(nam_dt, family) |>  
    summarise(gv = median(gv), pheno = median(pheno)) |> 
    dplyr::mutate(type = "Median") |> 
    select(family, gv, pheno, type)
  
  ###
  RP <- unique(rils_dt[[recur_parent_col]])
  DP <- unique(rils_dt[[donor_parent_col]])
  
  nfamilies = length(unique(nam_dt$family))
  source_dt <- if (is.null(parents_dt)) founders_dt else parents_dt
  
  parents_dt0 <- bind_rows(filter(source_dt, Founder_ID %in% DP),
                           filter(source_dt, Founder_ID == RP) |> 
                             dplyr::slice(rep(1, nfamilies)))
  parents_dt <- parents_dt0 |> 
    mutate(family = c(paste(RP, DP, sep = "X"), 
                      paste(RP, parents_dt0$Founder_ID[1:nfamilies], sep = "X")),
           gv = Trait1_GV, pheno = Trait1_Pheno,
           type = ifelse(Founder_ID == RP, "RP", "DP")) |> 
    select(family, gv, pheno, type)
  
  nam_dt_final <- rbind(nam_dt, nam_dt_summarized, parents_dt)
  
  #### Controlling plots parameters
  peak_pheno = 0
  peak_gv = 0
  for (f in unique(nam_dt_final$family)) {
    dens_pheno <- density(filter(nam_dt_final, type == "RIL" & family == f)$pheno, na.rm = T)
    dens_gv <- density(filter(nam_dt_final, type == "RIL"& family == f)$gv, na.rm = T)
    peak_pheno <- max(peak_pheno, max(dens_pheno$y))
    peak_gv <- max(peak_gv, max(dens_gv$y))
  }
  peak_pheno <<- peak_pheno + 0.05
  peak_gv <<- peak_gv + 0.05
  
  ### Now create plots 
  boris_theme <-theme_bw() + theme(axis.text = element_text(color = "black", size = 16),
                                   axis.title = element_text(colour = "black", face = "bold", size = 16),
                                   legend.text = element_text(colour = "black", size = 16),
                                   legend.position = "right", legend.direction = "vertical",
                                   legend.title = element_blank(),
                                   panel.grid = element_blank(),
                                   panel.background = element_rect(fill = "White", colour = NA),
                                   panel.border = element_rect(colour = "black"),
                                   text = element_text(colour = "black", size = 16))
  
  # pg_pheno <- ggplot()+
  #   geom_density(data = filter(nam_dt_final, type == "RIL"), aes(x = pheno), fill = "grey80") +
  #   geom_point(data =filter(nam_dt_final, type == "Median"), 
  #              aes(x = pheno, y = peak_pheno), size = 3, shape = 25, fill = "black", color  = "black")+
  #   geom_segment(data =filter(nam_dt_final, type == "RP"), 
  #                aes(x = pheno, xend = pheno, y = peak_pheno, yend = 0), linewidth = 1, color = "blue")+
  #   geom_segment(data =filter(nam_dt_final, type == "DP"), 
  #                aes(x = pheno, xend = pheno, y = peak_pheno, yend = 0), linewidth = 1, color = "darkred") +
  #   boris_theme + labs(x = paste0(trait_name, "(",trait_unit, ")"), y = "Density", title = "Phenotype") +
  #   scale_y_continuous(limits = c(0, peak_pheno), breaks = seq(0,peak_pheno,0.05)) +
  #   facet_wrap(~family, nrow = n_plot_rows, ncol = n_plot_cols)
  # 
  # pg_gv <- ggplot()+
  #   geom_density(data = filter(nam_dt_final, type == "RIL"), aes(x = gv), fill = "grey80") +
  #   geom_point(data =filter(nam_dt_final, type == "Median"), 
  #              aes(x = gv, y = peak_gv), size = 3, shape = 25, fill = "black", color  = "black")+
  #   geom_segment(data =filter(nam_dt_final, type == "RP"), 
  #                aes(x = gv, xend = gv, y = peak_gv, yend = 0), linewidth = 1, color = "blue")+
  #   geom_segment(data =filter(nam_dt_final, type == "DP"), 
  #                aes(x = gv, xend = gv, y = peak_gv, yend = 0), linewidth = 1, color = "darkred") +
  #   boris_theme + labs(x = paste0(trait_name, "(",trait_unit, ")"), y = "Density", title = "Genetic values") +
  #   scale_y_continuous(limits = c(0, peak_gv), breaks = seq(0,peak_gv,0.1)) +
  #   facet_wrap(~family, nrow = n_plot_rows, ncol = n_plot_cols)
  
  pg_pheno <- ggplot() +
    geom_density(
      data = filter(nam_dt_final, type == "RIL"),
      aes(x = pheno),
      fill = "grey80",
      color = NA
    ) +
    geom_point(
      data = filter(nam_dt_final, type == "Median"),
      aes(x = pheno, y = peak_pheno),
      size = 3, shape = 25, fill = "black", color = "black"
    ) +
    geom_segment(
      data = filter(nam_dt_final, type %in% c("RP","DP")),
      aes(x = pheno, xend = pheno, y = peak_pheno, yend = 0, color = type),
      linewidth = 1
    ) +
    boris_theme +
    labs(
      x = paste0(trait_name, " (", trait_unit, ")"),
      y = "Density",
      title = "Phenotype",
      color = "Founder Type"      # legend title
    ) +
    scale_color_manual(
      values = c(RP = "blue", DP = "darkred")
    ) +
    scale_y_continuous(
      limits = c(0, peak_pheno),
      breaks = seq(0, peak_pheno, 0.05)
    ) +
    facet_wrap(~family, nrow = n_plot_rows, ncol = n_plot_cols)
  
  pg_gv <- ggplot() +
    geom_density(
      data = filter(nam_dt_final, type == "RIL"),
      aes(x = gv),
      fill = "grey80",
      color = NA
    ) +
    geom_point(
      data = filter(nam_dt_final, type == "Median"),
      aes(x = gv, y = peak_gv),
      size = 3, shape = 25, fill = "black", color = "black"
    ) +
    geom_segment(
      data = filter(nam_dt_final, type %in% c("RP", "DP")),
      aes(x = gv, xend = gv, y = peak_gv, yend = 0, color = type),
      linewidth = 1
    ) +
    boris_theme +
    labs(
      x = paste0(trait_name, " (", trait_unit, ")"),
      y = "Density",
      title = "Genetic values",
      color = "Founder Type"       # Legend title
    ) +
    scale_color_manual(
      values = c(RP = "blue", DP = "darkred")
    ) +
    scale_y_continuous(
      limits = c(0, peak_gv),
      breaks = seq(0, peak_gv, 0.1)
    ) +
    facet_wrap(~family, nrow = n_plot_rows, ncol = n_plot_cols)
  
  ### Return a list of the plotted data and the plots.
  return(list(Final_data = nam_dt_final, plot_pheno = pg_pheno, plot_GV = pg_gv))
  
}