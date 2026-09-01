library(stats)

mat_sym_log <- read.table('PATH', sep = '\t', header = TRUE, check.names = FALSE) %>% as.matrix


dist_corr_log <- as.dist(1 - cor(t(mat_sym_log), use = "pairwise.complete.obs"))

mds_corr_log <- cmdscale(dist_corr_log, k = 2)


library(ggplot2)
library(ggrepel)

colnames(dist_corr_log) <- c('x', 'y')
mds_corr_log_df <- as.data.frame(mds_corr_log)
mds_corr_log_df$name <- rownames(mds_corr_log_df)

mds_corr_log_df$color <- sapply(mds_corr_log_df$name, function(x){
  if(grepl('\\(', x)){
    return(gsub('.*\\(', '', gsub('\\)', '', x)))
  } else return(x)
})

# color formatting
# mds_corr_log_df$color <- sapply(mds_corr_log_df$color, function(x){
#   result <- x
#   if(grepl('Central|Eastern|Middle|Northern|South|SouthEastern|Southern|Western', x)) result <- gsub('Central|Eastern|Middle|Northern|South|SouthEastern|Southern|Western', '', x)
#   if(result %in% c('America', 'Caribbean')) result <- 'Americas'
#   if(result %in% c('AustraliaandNewZealand', 'PacificIslands')) result <- 'Oceania'
#   return(result)
#   
# })


continent_colors <- c(
  "Africa"   = "#f58231",
  "Asia"     = "#e6194B",
  "Europe"   = "#3cb44b",
  "Americas" = "#000075",
  "Oceania"  = "#911eb4"
)


plt <- ggplot(as.data.frame(mds_corr_log_df), aes(x = x, y = y, label = name, colour = as.factor(color))) +
  geom_text_repel(
    size = 5,
    fontface = "bold",
    segment.color = "grey50",
    point.padding = 0.0,
    box.padding = 0.0,
    max.overlaps = 15,
    show.legend = FALSE  
  ) +
  geom_point(alpha = 0, size = 0) +
  labs(
    title = "Principal Coordinates Analysis of Relative Risks for ...",
    # subtitle = "",
    x = "Dimension 1",
    y = "Dimension 2",
    color = "Continental\nregion"
  ) +
  theme_minimal(base_size = 20) +
  theme(
    legend.title = element_text(size = 20), legend.text = element_text(size = 16),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
    legend.position = "right",
    panel.grid.major = element_line(colour = "grey90", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold")
  ) +  guides(colour = guide_legend(override.aes = list(shape = 15, size = 10, alpha = 1))) +
  scale_colour_manual(values = continent_colors)



ggsave('PATH TO SAVE', plot = plt, width = 12, height = 15, dpi = 200, bg = "white")


write.table(mds_corr_log_df, 'PATH TO SAVE 2', sep = '\t', row.names = FALSE, quote = FALSE)











