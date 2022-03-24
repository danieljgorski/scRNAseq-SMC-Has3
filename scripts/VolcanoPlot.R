VolcanoPlot <- function (deg_data_frame, title) {

  # Calculate -log10(p-adj)
  deg <- deg_data_frame
  deg$neglog10p <- -(log10(deg$p_val_adj))
  
  # filter significant
  deg_up <- deg[deg$regulation == "Up",]
  deg_up_sig <- deg_up[deg_up$p_val_adj < 0.01,]
  deg_up_sig_fc <- deg_up_sig[deg_up_sig$avg_logFC > 0.25,]
  deg_up_label <- deg_up_sig_fc %>% slice_max(neglog10p, n = 20)
  deg_down <- deg[deg$regulation == "Down",]
  deg_down_sig <- deg_down[deg_down$p_val_adj < 0.01,]
  deg_down_sig_fc <- deg_down_sig[deg_down_sig$avg_logFC < -0.25,]
  deg_down_label <- deg_down_sig_fc %>% slice_max(neglog10p, n = 20)
  
  # volcano plot
  vp <- ggplot(deg,
               aes(x = avg_logFC, 
                   y = neglog10p)) 
  vc <- vp +
    geom_point(colour = "grey") +
    geom_point(data = deg_up_sig, colour="#eb2a0e") +
    geom_point(data = deg_down_sig, colour="#2664ad") +
    geom_text_repel(data = deg_up_label, 
                    aes(label=gene), 
                    colour = "#eb2a0e", 
                    fontface = "italic", 
                    size = 3,
                    force = 2,
                    nudge_x  = max(abs(deg$avg_logFC)),
                    nudge_y = median(abs(deg$neglog10p)),
                    direction = "both",
                    segment.size = 0.1,
                    segment.alpha = 0.3,
                    segment.curvature =0,
                    max.overlaps = Inf) +
    geom_text_repel(data = deg_down_label, 
                    aes(label=gene), 
                    colour = "#2664ad", 
                    fontface = "italic", 
                    size = 3,
                    force = 2,
                    nudge_x  = -(max(abs(deg$avg_logFC))),
                    nudge_y = median(abs(deg$neglog10p)),
                    direction = "both",
                    segment.size = 0.1,
                    segment.alpha = 0.3,
                    segment.curvature = 0,
                    max.overlaps = Inf) +
    xlab("avg logFC") +
    ylab("-log10 (p-adj)") +
    xlim(-max(abs(deg$avg_logFC)), 
         max(abs(deg$avg_logFC))) +
    geom_hline(yintercept = -log10(0.01), 
               linetype="dashed", 
               color = "grey", 
               size = 0.1) +
    geom_vline(xintercept = -0.25, 
               linetype="dashed", 
               color = "grey", 
               size = 0.1) +
    geom_vline(xintercept = 0.25, 
               linetype="dashed", 
               color = "grey",
               size = 0.1) +
    theme_bw() + 
    theme(panel.grid.major =  element_blank(),
          panel.grid.minor = element_blank(),
          axis.title = element_text(size = 14)) +
    ggtitle(paste0(title))
  print(vc)
}


