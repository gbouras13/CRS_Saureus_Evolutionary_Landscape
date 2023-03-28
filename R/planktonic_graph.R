library(tidyverse)

# Load data
df_prop <- read.csv("df_prop.csv")

# Set factor levels for Susceptibility
df_prop$Susceptibility <- factor(df_prop$Susceptibility, levels = c("susceptible", "intermediate", "resistant"))

# Create arc data frame
r <- 0.15
t <- seq(0, 179, by = 1) * pi / 180
x <- r * cos(t)
y <- r * 5 * sin(t)
arc_df <- data.frame(Group = x + 1.5, Value = y + 102)

# Function to create plot
create_plot <- function(data, title) {
  plot <- ggplot(data, aes(timepoint, proportion)) +
    geom_col(aes(fill = Susceptibility, colour = Susceptibility)) +
    geom_text(x = 1.5, y = 105, label = "N.S") +
    facet_grid(~ antibiotic) +
    geom_line(data = arc_df, aes(Group, Value)) +
    theme_classic(base_size = 16) +
    # theme(strip.placement = "outside") +
    scale_fill_manual(values = c("#FDDBC7FF","#D6604DFF","#B2182BFF")) +
    scale_color_manual(values = c("#FDDBC7FF","#D6604DFF","#B2182BFF")) +
    ylab("Susceptibility (%)") +
    xlab("") +
    ggtitle(title)
  return(plot)
}

# Create plots
diff_strain_plot <- df_prop %>%
  filter(same_strain == 0) %>%
  create_plot(title = "Different strain pairs")

same_strain_plot <- df_prop %>%
  filter(same_strain == 1) %>%
  create_plot(title = "Same strain pairs")

# Arrange plots
combined_plot <- ggpubr::ggarrange(diff_strain_plot, same_strain_plot, ncol = 1, common.legend = TRUE, legend = "top")

# Display plot
combined_plot



svglite:: svglite(filename =paste0("Figures/","planktonic.svg"), width = 8, height = 8)
combined_plot
dev.off()


png(filename =paste0("Figures/", "planktonic.png"), units = "in", width = 5.77*2, height = 10, res = 1200,  type='cairo') 
combined_plot
dev.off()


