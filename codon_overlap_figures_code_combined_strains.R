library(ggplot2)
library(readr)
library(dplyr)
library(multcomp)
library(ggpubr)
library(ggsignif)
library(car)
library(tidyr)
library(MASS)
library(gridExtra)
library(gtable)
library(grid)

# Reading datasets from CSV files into data frames and adding a source column
wglo <- read_csv('D:/Users/jkm78/Spyder/test_df_87_overlap_prediction_from_uniform_knowns.csv') %>% mutate(source = 'W. glossinidia (23.6 GC)')
sflo <- read_csv('D:/Users/jkm78/Spyder/test_df_88_overlap_prediction_from_uniform_knowns.csv') %>% mutate(source = 'S. floricola (24.8 GC)')
aman <- read_csv('D:/Users/jkm78/Spyder/test_df_89_overlap_prediction_from_uniform_knowns.csv') %>% mutate(source = 'A. manzaensis (31.4 GC)')
psal <- read_csv('D:/Users/jkm78/Spyder/test_df_90_overlap_prediction_from_uniform_knowns.csv') %>% mutate(source = 'P. saltans (37.4 GC)')
cmaq <- read_csv('D:/Users/jkm78/Spyder/test_df_91_overlap_prediction_from_uniform_knowns.csv') %>% mutate(source = 'C. maquilingensis (43.9 GC)')
sfle <- read_csv('D:/Users/jkm78/Spyder/test_df_92_overlap_prediction_from_uniform_knowns.csv') %>% mutate(source = 'S. flexneri (50.9 GC)')
lfer <- read_csv('D:/Users/jkm78/Spyder/test_df_93_overlap_prediction_from_uniform_knowns.csv') %>% mutate(source = 'L. fermentum (53.0 GC)')
gkil <- read_csv('D:/Users/jkm78/Spyder/test_df_94_overlap_prediction_from_uniform_knowns.csv') %>% mutate(source = 'G. kilaueensis (61.4 GC)')
obac <- read_csv('D:/Users/jkm78/Spyder/test_df_95_overlap_prediction_from_uniform_knowns.csv') %>% mutate(source = 'O. bacterium (64.5 GC)')
rdep <- read_csv('D:/Users/jkm78/Spyder/test_df_96_overlap_prediction_from_uniform_knowns.csv') %>% mutate(source = 'R. depolymerans (67.3 GC)')
afri <- read_csv('D:/Users/jkm78/Spyder/test_df_97_overlap_prediction_from_uniform_knowns.csv') %>% mutate(source = 'A. friuliensis (70.8 GC)')
gobs <- read_csv('D:/Users/jkm78/Spyder/test_df_98_overlap_prediction_from_uniform_knowns.csv') %>% mutate(source = 'G. obscurus (74.1 GC)')


# Combining all datasets into one dataframe
combined_df <- bind_rows(wglo,
                         sflo,
                         aman,
                         psal,
                         cmaq,
                         sfle,
                         lfer,
                         gkil,
                         obac,
                         rdep,
                         afri,
                         gobs)

# Filter the combined dataframe based on your conditions
filtered_df <- combined_df %>% 
  filter(length_2 > 4, length_2 == length_actual, converted_overlap != "none")

# Filter the combined dataframe based on the specified length ranges
filtered_df <- combined_df %>%
  filter((length_2 >= 6 & length_2 <= 12) | 
           (length_2 >= 52 & length_2 <= 54) | 
           (length_2 >= 97 & length_2 <= 99),
         length_2 == length_actual, 
         converted_overlap != "none")

# Prepare summary data for bar chart
summary <- filtered_df %>% 
  group_by(length_actual, source) %>%
  summarise(count = n(), .groups = 'drop') %>%
  mutate(percentage_from_total = (count / sum(count)) * 100)

table(new_df$length_actual)

# Ensure the order of the groups in the plot
ordered_sources <- c('W. glossinidia (23.6 GC)', 
                     'S. floricola (24.8 GC)', 
                     'A. manzaensis (31.4 GC)', 
                     'P. saltans (37.4 GC)', 
                     'C. maquilingensis (43.9 GC)', 
                     'S. flexneri (50.9 GC)', 
                     'L. fermentum (53.0 GC)', 
                     'G. kilaueensis (61.4 GC)',
                     'O. bacterium (64.5 GC)',
                     'R. depolymerans (67.3 GC)',
                     'A. friuliensis (70.8 GC)',
                     'G. obscurus (74.1 GC)')

filtered_df$source <- factor(filtered_df$source, levels = ordered_sources)

# Plotting
p <- ggplot() +
  geom_boxplot(data = filtered_df, aes(x=factor(length_2), y=normalized_alignment_score_2, fill=source)) +
  labs(
    title = NULL,
    x = "Overlap Length",
    y = "Alignment Score",
    fill = NULL
  ) +
  
  # Primary y-axis
  scale_y_continuous(
    limits = c(0, 1)
  ) +
  
  theme(
    axis.text.x = element_text(size = 10,angle = 45, hjust = 1.0, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 10, color = "black"),
    plot.title = element_text(hjust = 0.5),
    legend.text = element_text(size = 10, color = "black"),  # Adjusted legend text size
    legend.title = element_text(size = 16, color = "black"),
    legend.position = "right",
    text = element_text(color = "black")
  )

print(p)

######### Generate a combined figure with one on top, and two on the bottom

# Filter the combined dataframe based on the specified length ranges
filtered_df_top <- combined_df %>%
  filter((length_2 >= 6 & length_2 <= 18),
         length_2 == length_actual, 
         converted_overlap != "none")

# Ensure the order of the groups in the plot
ordered_sources <- c('W. glossinidia (23.6 GC)', 
                     'S. floricola (24.8 GC)', 
                     'A. manzaensis (31.4 GC)', 
                     'P. saltans (37.4 GC)', 
                     'C. maquilingensis (43.9 GC)', 
                     'S. flexneri (50.9 GC)', 
                     'L. fermentum (53.0 GC)', 
                     'G. kilaueensis (61.4 GC)',
                     'O. bacterium (64.5 GC)',
                     'R. depolymerans (67.3 GC)',
                     'A. friuliensis (70.8 GC)',
                     'G. obscurus (74.1 GC)')

filtered_df_top$source <- factor(filtered_df_top$source, levels = ordered_sources)

# Filtered data for the bottom plots
filtered_df_bottom_left <- combined_df %>%
  filter((length_2 >= 52 & length_2 <= 54),
         length_2 == length_actual, 
         converted_overlap != "none")

filtered_df_bottom_right <- combined_df %>%
  filter((length_2 >= 97 & length_2 <= 99),
         length_2 == length_actual, 
         converted_overlap != "none")

# Ensure the same factor levels for "source" across all filtered datasets
filtered_df_bottom_left$source <- factor(filtered_df_bottom_left$source, levels = ordered_sources)
filtered_df_bottom_right$source <- factor(filtered_df_bottom_right$source, levels = ordered_sources)

# Plotting the top plot
top_plot <- ggplot(data = filtered_df_top, aes(x = factor(length_2), y = normalized_alignment_score_2, fill = source)) +
  geom_boxplot() +
  labs(title = NULL, x = "Overlap Length", y = "Alignment Score", fill = NULL) +
  scale_y_continuous(limits = c(0, 1)) +
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1.0, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 14, color = "black"),
    axis.title.y = element_text(size = 14, color = "black"),
    plot.title = element_text(hjust = 0.5),
    legend.text = element_text(size = 10, color = "black"),  # Adjusted legend text size
    legend.title = element_text(size = 16, color = "black"),
    text = element_text(color = "black")
  )

# Plotting the bottom left plot
bottom_left_plot <- ggplot(data = filtered_df_bottom_left, aes(x = factor(length_2), y = normalized_alignment_score_2, fill = source)) +
  geom_boxplot() +
  labs(title = NULL, x = "Overlap Length", y = "Alignment Score", fill = NULL) +
  scale_y_continuous(limits = c(0, 1)) +
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1.0, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 14, color = "black"),
    axis.title.y = element_text(size = 14, color = "black"),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    text = element_text(color = "black")
  )

# Plotting the bottom right plot
bottom_right_plot <- ggplot(data = filtered_df_bottom_right, aes(x = factor(length_2), y = normalized_alignment_score_2, fill = source)) +
  geom_boxplot() +
  labs(title = NULL, x = "Overlap Length", y = NULL, fill = NULL) +
  scale_y_continuous(limits = c(0, 1)) +
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1.0, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 14, color = "black"),
    axis.title.y = element_text(size = 14, color = "black"),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    text = element_text(color = "black")
  )

# Adjust the top plot to create space for the legend at the bottom but don't display it
top_plot <- top_plot + theme(legend.position = "none")

# Combine the bottom plots into a single grob
bottom_plots <- arrangeGrob(bottom_left_plot, bottom_right_plot, ncol = 2)

# Extract and adjust the legend to have two rows
legend_plot <- top_plot + theme(legend.position = "bottom") + guides(fill = guide_legend(nrow = 4, byrow = TRUE))
legend <- gtable::gtable_filter(ggplotGrob(legend_plot), "guide-box")

spacing_factor = 0.25  # Adjust this factor to increase or decrease the space
adjusted_heights_for_plots = c(1, 2 + spacing_factor)

# Now, arrange the plots and legend with the adjusted spacing
final_plot <- arrangeGrob(
  top_plot,
  bottom_plots,
  legend,
  heights = unit(c(adjusted_heights_for_plots, 0.5), "null")  # Adjusted heights for spacing
)

# Draw the plot
grid.draw(final_plot) #save to pdf as 7.5 x 8.5 in

# Since we now ignore top_plot, extract the legend directly from bottom_left_plot or bottom_right_plot (assuming they have the same legend)
# Let's use bottom_left_plot for extracting legend. First, create a version of bottom_left_plot that includes the legend at the bottom.
legend_plot <- bottom_left_plot + theme(legend.position = "bottom") + guides(fill = guide_legend(nrow = 3, byrow = TRUE))
legend <- gtable::gtable_filter(ggplotGrob(legend_plot), "guide-box")

# Combine bottom_left_plot and bottom_right_plot without any legends displayed within them
bottom_left_plot <- bottom_left_plot + theme(legend.position = "none")
bottom_right_plot <- bottom_right_plot + theme(legend.position = "none")
bottom_plots <- arrangeGrob(bottom_left_plot, bottom_right_plot, ncol = 2)

# Arrange the bottom plots and the extracted legend
# Adjust the heights to allocate appropriate space to plots and legend
final_plot <- arrangeGrob(
  bottom_plots,
  legend,
  heights = unit(c(1, 0.2), "null")  # Height adjustments: larger space for plots, less for legend
)

# Draw the plot
grid.draw(final_plot)

######## stats tests

################

# Assuming 'filtered_df' is your dataset and 'normalized_alignment_score_2' is the variable of interest
length_97_df <- filtered_df %>% filter(length_2 == 97)

anova_result <- aov(normalized_alignment_score_2 ~ source, data = length_97_df)
summary(anova_result)

# Faceted Histogram
ggplot(length_97_df, aes(x = normalized_alignment_score_2)) + 
  geom_histogram(bins = 30, fill = "blue", alpha = 0.7) + 
  facet_wrap(~ source) +
  theme_minimal() +
  labs(title = "Histogram of Normalized Alignment Score by Source", 
       x = "Normalized Alignment Score", 
       y = "Count")

#Note, these data do not all appear normally distributed

# Q-Q plots for each group
ggqqplot(length_97_df, "normalized_alignment_score_2", facet.by = "source", 
         ggtheme = theme_bw(), color = "source")

# Shapiro-Wilk test for each group
grouped_normality_test <- length_97_df %>%
  group_by(source) %>%
  do(shapiro_test = shapiro.test(.$normalized_alignment_score_2))

print(grouped_normality_test)

grouped_normality_test$shapiro_test

#yep, shapiro test indicates non-normality

# Levene test to check for homogeneity of variance
levene_test <- leveneTest(normalized_alignment_score_2 ~ source, data = length_97_df)
print(levene_test)

#levene test indicates non-homogenous variance

#Using non-parametric Kruskal-Wallis test, followed by pairwise Wilcoxon rank sum tests (also known as the Mann-Whitney U test),
#with bonferroni correction applied to adjust p-values (as this is a conservative test)

kruskal_result <- kruskal.test(normalized_alignment_score_2 ~ source, data = length_97_df)
print(kruskal_result)

pairwise_results <- pairwise.wilcox.test(length_97_df$normalized_alignment_score_2, length_97_df$source, 
                     p.adjust.method = "bonferroni")

print(pairwise_results)

################

# Assuming 'filtered_df' is your dataset and 'normalized_alignment_score_2' is the variable of interest
length_98_df <- filtered_df %>% filter(length_2 == 98)

# Faceted Histogram
ggplot(length_98_df, aes(x = normalized_alignment_score_2)) + 
  geom_histogram(bins = 30, fill = "blue", alpha = 0.7) + 
  facet_wrap(~ source) +
  theme_minimal() +
  labs(title = "Histogram of Normalized Alignment Score by Source", 
       x = "Normalized Alignment Score", 
       y = "Count")

# Q-Q plots for each group
ggqqplot(length_98_df, "normalized_alignment_score_2", facet.by = "source", 
         ggtheme = theme_bw(), color = "source")

# Shapiro-Wilk test for each group
grouped_normality_test <- length_98_df %>%
  group_by(source) %>%
  do(shapiro_test = shapiro.test(.$normalized_alignment_score_2))

print(grouped_normality_test)

grouped_normality_test$shapiro_test

#yep, shapiro test indicates non-normality

# Levene test to check for homogeneity of variance
levene_test <- leveneTest(normalized_alignment_score_2 ~ source, data = length_98_df)
print(levene_test)

#levene test indicates non-homogenous variance

#Using non-parametric Kruskal-Wallis test, followed by pairwise Wilcoxon rank sum tests (also known as the Mann-Whitney U test),
#with bonferroni correction applied to adjust p-values (as this is a conservative test)

kruskal_result <- kruskal.test(normalized_alignment_score_2 ~ source, data = length_98_df)
print(kruskal_result)

pairwise_results <- pairwise.wilcox.test(length_98_df$normalized_alignment_score_2, length_98_df$source, 
                                         p.adjust.method = "bonferroni")

print(pairwise_results)

################

# Assuming 'filtered_df' is your dataset and 'normalized_alignment_score_2' is the variable of interest
length_99_df <- filtered_df %>% filter(length_2 == 99)

# Faceted Histogram
ggplot(length_99_df, aes(x = normalized_alignment_score_2)) + 
  geom_histogram(bins = 30, fill = "blue", alpha = 0.7) + 
  facet_wrap(~ source) +
  theme_minimal() +
  labs(title = "Histogram of Normalized Alignment Score by Source", 
       x = "Normalized Alignment Score", 
       y = "Count")

# Q-Q plots for each group
ggqqplot(length_98_df, "normalized_alignment_score_2", facet.by = "source", 
         ggtheme = theme_bw(), color = "source")

# Shapiro-Wilk test for each group
grouped_normality_test <- length_99_df %>%
  group_by(source) %>%
  do(shapiro_test = shapiro.test(.$normalized_alignment_score_2))

print(grouped_normality_test)

grouped_normality_test$shapiro_test

#yep, shapiro test indicates non-normality

# Levene test to check for homogeneity of variance
levene_test <- leveneTest(normalized_alignment_score_2 ~ source, data = length_99_df)
print(levene_test)

#levene test indicates non-homogenous variance

#Using non-parametric Kruskal-Wallis test, followed by pairwise Wilcoxon rank sum tests (also known as the Mann-Whitney U test),
#with bonferroni correction applied to adjust p-values (as this is a conservative test)

kruskal_result <- kruskal.test(normalized_alignment_score_2 ~ source, data = length_99_df)
print(kruskal_result)

pairwise_results <- pairwise.wilcox.test(length_99_df$normalized_alignment_score_2, length_99_df$source, 
                                         p.adjust.method = "bonferroni")

print(pairwise_results)

# Assuming 'filtered_df' is your dataset and 'normalized_alignment_score_2' is the variable of interest
length_97_df <- filtered_df %>% filter(length_2 == 97)

anova_result <- aov(normalized_alignment_score_2 ~ source, data = length_97_df)
summary(anova_result)

# Faceted Histogram
ggplot(length_97_df, aes(x = normalized_alignment_score_2)) + 
  geom_histogram(bins = 30, fill = "blue", alpha = 0.7) + 
  facet_wrap(~ source) +
  theme_minimal() +
  labs(title = "Histogram of Normalized Alignment Score by Source", 
       x = "Normalized Alignment Score", 
       y = "Count")

#Note, these data do not all appear normally distributed

# Q-Q plots for each group
ggqqplot(length_97_df, "normalized_alignment_score_2", facet.by = "source", 
         ggtheme = theme_bw(), color = "source")

# Shapiro-Wilk test for each group
grouped_normality_test <- length_97_df %>%
  group_by(source) %>%
  do(shapiro_test = shapiro.test(.$normalized_alignment_score_2))

print(grouped_normality_test)

grouped_normality_test$shapiro_test

#yep, shapiro test indicates non-normality

# Levene test to check for homogeneity of variance
levene_test <- leveneTest(normalized_alignment_score_2 ~ source, data = length_97_df)
print(levene_test)

#levene test indicates non-homogenous variance

#Using non-parametric Kruskal-Wallis test, followed by pairwise Wilcoxon rank sum tests (also known as the Mann-Whitney U test),
#with bonferroni correction applied to adjust p-values (as this is a conservative test)

kruskal_result <- kruskal.test(normalized_alignment_score_2 ~ source, data = length_97_df)
print(kruskal_result)

pairwise_results <- pairwise.wilcox.test(length_97_df$normalized_alignment_score_2, length_97_df$source, 
                     p.adjust.method = "bonferroni")

print(pairwise_results)

################

# Assuming 'filtered_df' is your dataset and 'normalized_alignment_score_2' is the variable of interest
length_98_df <- filtered_df %>% filter(length_2 == 98)

# Faceted Histogram
ggplot(length_98_df, aes(x = normalized_alignment_score_2)) + 
  geom_histogram(bins = 30, fill = "blue", alpha = 0.7) + 
  facet_wrap(~ source) +
  theme_minimal() +
  labs(title = "Histogram of Normalized Alignment Score by Source", 
       x = "Normalized Alignment Score", 
       y = "Count")

# Q-Q plots for each group
ggqqplot(length_98_df, "normalized_alignment_score_2", facet.by = "source", 
         ggtheme = theme_bw(), color = "source")

# Shapiro-Wilk test for each group
grouped_normality_test <- length_98_df %>%
  group_by(source) %>%
  do(shapiro_test = shapiro.test(.$normalized_alignment_score_2))

print(grouped_normality_test)

grouped_normality_test$shapiro_test

#yep, shapiro test indicates non-normality

# Levene test to check for homogeneity of variance
levene_test <- leveneTest(normalized_alignment_score_2 ~ source, data = length_98_df)
print(levene_test)

#levene test indicates non-homogenous variance

#Using non-parametric Kruskal-Wallis test, followed by pairwise Wilcoxon rank sum tests (also known as the Mann-Whitney U test),
#with bonferroni correction applied to adjust p-values (as this is a conservative test)

kruskal_result <- kruskal.test(normalized_alignment_score_2 ~ source, data = length_98_df)
print(kruskal_result)

pairwise_results <- pairwise.wilcox.test(length_98_df$normalized_alignment_score_2, length_98_df$source, 
                                         p.adjust.method = "bonferroni")

print(pairwise_results)

################

# Assuming 'filtered_df' is your dataset and 'normalized_alignment_score_2' is the variable of interest
length_99_df <- filtered_df %>% filter(length_2 == 99)

# Faceted Histogram
ggplot(length_99_df, aes(x = normalized_alignment_score_2)) + 
  geom_histogram(bins = 30, fill = "blue", alpha = 0.7) + 
  facet_wrap(~ source) +
  theme_minimal() +
  labs(title = "Histogram of Normalized Alignment Score by Source", 
       x = "Normalized Alignment Score", 
       y = "Count")

# Q-Q plots for each group
ggqqplot(length_98_df, "normalized_alignment_score_2", facet.by = "source", 
         ggtheme = theme_bw(), color = "source")

# Shapiro-Wilk test for each group
grouped_normality_test <- length_99_df %>%
  group_by(source) %>%
  do(shapiro_test = shapiro.test(.$normalized_alignment_score_2))

print(grouped_normality_test)

grouped_normality_test$shapiro_test

#yep, shapiro test indicates non-normality

# Levene test to check for homogeneity of variance
levene_test <- leveneTest(normalized_alignment_score_2 ~ source, data = length_99_df)
print(levene_test)

#levene test indicates non-homogenous variance

#Using non-parametric Kruskal-Wallis test, followed by pairwise Wilcoxon rank sum tests (also known as the Mann-Whitney U test),
#with bonferroni correction applied to adjust p-values (as this is a conservative test)

kruskal_result <- kruskal.test(normalized_alignment_score_2 ~ source, data = length_99_df)
print(kruskal_result)

pairwise_results <- pairwise.wilcox.test(length_99_df$normalized_alignment_score_2, length_99_df$source, 
                                         p.adjust.method = "bonferroni")

print(pairwise_results)

#############

# Assuming 'filtered_df' is your dataset and 'normalized_alignment_score_2' is the variable of interest
length_52_df <- filtered_df %>% filter(length_2 == 52)

# Faceted Histogram
ggplot(length_52_df, aes(x = normalized_alignment_score_2)) + 
  geom_histogram(bins = 30, fill = "blue", alpha = 0.7) + 
  facet_wrap(~ source) +
  theme_minimal() +
  labs(title = "Histogram of Normalized Alignment Score by Source for Length 52", 
       x = "Normalized Alignment Score", 
       y = "Count")

# Q-Q plots for each group, correcting the reference dataframe
ggqqplot(length_52_df, "normalized_alignment_score_2", facet.by = "source", 
         ggtheme = theme_bw(), color = "source")

# Shapiro-Wilk test for each group
grouped_normality_test <- length_52_df %>%
  group_by(source) %>%
  summarise(shapiro_test_p_value = shapiro.test(normalized_alignment_score_2)$p.value)

print(grouped_normality_test)

# Levene test to check for homogeneity of variance
levene_test <- leveneTest(normalized_alignment_score_2 ~ source, data = length_52_df)
print(levene_test)

# Using non-parametric Kruskal-Wallis test, followed by pairwise Wilcoxon rank sum tests
kruskal_result <- kruskal.test(normalized_alignment_score_2 ~ source, data = length_52_df)
print(kruskal_result)

# Performing pairwise comparisons using different methods for p-value adjustment


pairwise_results_bonferroni <- pairwise.wilcox.test(length_52_df$normalized_alignment_score_2, length_52_df$source, 
                                                    p.adjust.method = "bonferroni")
print(pairwise_results_bonferroni)

#################

# Assuming 'filtered_df' is your dataset and 'normalized_alignment_score_2' is the variable of interest
length_53_df <- filtered_df %>% filter(length_2 == 53)

# Faceted Histogram
ggplot(length_53_df, aes(x = normalized_alignment_score_2)) + 
  geom_histogram(bins = 30, fill = "blue", alpha = 0.7) + 
  facet_wrap(~ source) +
  theme_minimal() +
  labs(title = "Histogram of Normalized Alignment Score by Source for Length 53", 
       x = "Normalized Alignment Score", 
       y = "Count")

# Q-Q plots for each group
ggqqplot(length_53_df, "normalized_alignment_score_2", facet.by = "source", 
         ggtheme = theme_bw(), color = "source")

# Shapiro-Wilk test for each group
grouped_normality_test <- length_53_df %>%
  group_by(source) %>%
  summarise(shapiro_test_p_value = shapiro.test(normalized_alignment_score_2)$p.value)

print(grouped_normality_test)

# Levene test to check for homogeneity of variance
levene_test <- leveneTest(normalized_alignment_score_2 ~ source, data = length_53_df)
print(levene_test)

# Kruskal-Wallis test
kruskal_result <- kruskal.test(normalized_alignment_score_2 ~ source, data = length_53_df)
print(kruskal_result)

# Pairwise Wilcoxon tests with multiple p-value adjustments
pairwise_results_bonferroni <- pairwise.wilcox.test(length_53_df$normalized_alignment_score_2, length_53_df$source, 
                                                    p.adjust.method = "bonferroni")
print(pairwise_results_bonferroni)


###########

# Assuming 'filtered_df' is your dataset and 'normalized_alignment_score_2' is the variable of interest
length_54_df <- filtered_df %>% filter(length_2 == 54)

# Faceted Histogram
ggplot(length_54_df, aes(x = normalized_alignment_score_2)) + 
  geom_histogram(bins = 30, fill = "blue", alpha = 0.7) + 
  facet_wrap(~ source) +
  theme_minimal() +
  labs(title = "Histogram of Normalized Alignment Score by Source for Length 54", 
       x = "Normalized Alignment Score", 
       y = "Count")

# Q-Q plots for each group
ggqqplot(length_54_df, "normalized_alignment_score_2", facet.by = "source", 
         ggtheme = theme_bw(), color = "source")

# Shapiro-Wilk test for each group
grouped_normality_test <- length_54_df %>%
  group_by(source) %>%
  summarise(shapiro_test_p_value = shapiro.test(normalized_alignment_score_2)$p.value)

print(grouped_normality_test)

# Levene test to check for homogeneity of variance
levene_test <- leveneTest(normalized_alignment_score_2 ~ source, data = length_54_df)
print(levene_test)

# Kruskal-Wallis test
kruskal_result <- kruskal.test(normalized_alignment_score_2 ~ source, data = length_54_df)
print(kruskal_result)

# Pairwise Wilcoxon tests with multiple p-value adjustments
pairwise_results_bonferroni <- pairwise.wilcox.test(length_54_df$normalized_alignment_score_2, length_54_df$source, 
                                                    p.adjust.method = "bonferroni")
print(pairwise_results_bonferroni)


