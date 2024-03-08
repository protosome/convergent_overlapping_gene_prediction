library(ggplot2)
library(readr)
library(gridExtra)
library(dplyr)
library(e1071)


# Read each CSV file into its own dataframe
df_87 <- read_csv('D:/Users/jkm78/Spyder/test_df_87_overlap_prediction_from_uniform_knowns.csv')
df_88 <- read_csv('D:/Users/jkm78/Spyder/test_df_88_overlap_prediction_from_uniform_knowns.csv')
df_89 <- read_csv('D:/Users/jkm78/Spyder/test_df_89_overlap_prediction_from_uniform_knowns.csv')
df_90 <- read_csv('D:/Users/jkm78/Spyder/test_df_90_overlap_prediction_from_uniform_knowns.csv')
df_91 <- read_csv('D:/Users/jkm78/Spyder/test_df_91_overlap_prediction_from_uniform_knowns.csv')
df_92 <- read_csv('D:/Users/jkm78/Spyder/test_df_92_overlap_prediction_from_uniform_knowns.csv')
df_93 <- read_csv('D:/Users/jkm78/Spyder/test_df_93_overlap_prediction_from_uniform_knowns.csv')
df_94 <- read_csv('D:/Users/jkm78/Spyder/test_df_94_overlap_prediction_from_uniform_knowns.csv')
df_95 <- read_csv('D:/Users/jkm78/Spyder/test_df_95_overlap_prediction_from_uniform_knowns.csv')
df_96 <- read_csv('D:/Users/jkm78/Spyder/test_df_96_overlap_prediction_from_uniform_knowns.csv')
df_97 <- read_csv('D:/Users/jkm78/Spyder/test_df_97_overlap_prediction_from_uniform_knowns.csv')
df_98 <- read_csv('D:/Users/jkm78/Spyder/test_df_98_overlap_prediction_from_uniform_knowns.csv')

# Combine all dataframes into one using rbind
new_df <- rbind(df_87, df_88, df_89, df_90, df_91, df_92, 
                df_93, df_94, df_95, df_96, df_97, df_98)


##################
# Analyzing the number of correct predictions (overlap or no overlap)

# Load the data
#df <- read.csv("D:/Users/jkm78/Spyder/test_df_36_overlap_prediction_from_uniform_knowns.csv")
df <- new_df

# Calculate the total number of rows where length_actual is 1
total_count <- nrow(subset(new_df, length_actual == 6)) #selecting length_actual representing group counts
total_count

# 3. Filter the data based on your conditions
filtered_df <- df %>% filter(length_actual > 4 & length_2 == 1)

# 4. Calculate the counts and percentages
summary <- filtered_df %>% 
  group_by(length_actual) %>%
  summarise(count = n()) %>%
  mutate(percentage_from_total = (count / total_count) * 100) %>%
  left_join(df %>% group_by(length_actual) %>% summarise(total_count = n()), by = "length_actual")

# 5. Display the table
print(summary, n = 29)

summary(summary$percentage_from_total)

#########
#combined bar and box plots
#########


# Base plot
p <- ggplot() +
  # Boxplot
  geom_boxplot(data = filtered_length_4, aes(x=factor(length_2), y=normalized_alignment_score_2, fill=factor(length_2))) +
  
  # Bar chart
  geom_col(data=summary, aes(x=factor(length_actual), y=percentage_from_total/100, fill=factor(length_actual)), position="dodge", width=0.5) +
  
  labs(
    #title = "Model Performance for Different Overlap Sequence Lengths",
    x = "Overlap Length",
    y = "Alignment Score",
    fill = "test"
  ) +
  
  # Primary y-axis
  scale_y_continuous(
    limits = c(0, 1),
    sec.axis = sec_axis(~.*100, name = "Percent Without Identification")
  ) +
  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1.0, size = 18),
    axis.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22, vjust = 2),
    axis.title.y.right = element_text(vjust = 1.7),
    plot.title = element_text(size = 18, hjust = 0.5),
    legend.text = element_text(18),
    legend.title = element_text(size = 18),
    legend.position = "none"
  )

# Modify x-axis labels to show at intervals of 2:
p + scale_x_discrete(breaks = unique(as.character(seq(from = min(as.numeric(as.character(filtered_length_4$length_2))), 
                                                      to = max(as.numeric(as.character(filtered_length_4$length_2))), by = 4))))



############################################################
### Summarizing the data by length_2, calculating the mean, median, and number of observations for each length_2
############################################################

summary_by_length <- filtered_length_4 %>%
  group_by(length_2) %>%
  summarise(
    mean_score = mean(normalized_alignment_score_2, na.rm = TRUE),
    median_score = median(normalized_alignment_score_2, na.rm = TRUE),
    sd_score = sd(normalized_alignment_score_2, na.rm = TRUE),
    n = n(),
    sem = sd_score / sqrt(n), # Standard Error of the Mean
    ci_lower = mean_score - (1.96 * sem), # Lower bound of the 95% CI
    ci_upper = mean_score + (1.96 * sem), # Upper bound of the 95% CI
    min_score = min(normalized_alignment_score_2, na.rm = TRUE), # Minimum score
    max_score = max(normalized_alignment_score_2, na.rm = TRUE), # Maximum score
    iqr_score = IQR(normalized_alignment_score_2, na.rm = TRUE), # Interquartile Range
    skewness = e1071::skewness(normalized_alignment_score_2, na.rm = TRUE), # Skewness
    kurtosis = e1071::kurtosis(normalized_alignment_score_2, na.rm = TRUE) # Kurtosis
  )

# View the summarized data
print(summary_by_length)

#build a function that will add phase information based on length
assign_phase <- function(value) {
  remainder <- value %% 3
  if (remainder == 0 || value == 6) {
    return(0)
  } else if (remainder == 1 || value == 4) {
    return(2)
  } else {
    return(1)
  }
}

summary_by_length$phase <- as.integer(lapply(summary_by_length$length_2, assign_phase))

####
#### Binning the mean_score data to summarize
####

# Function to create bins of length 10 with the last bin having at least 5 entries
# Function to create custom bins
create_bins <- function(lengths) {
  bin_breaks <- c(0, 18, 40, 66, Inf)  # Define custom breaks
  bins <- cut(lengths, breaks = bin_breaks, labels = FALSE)
  return(bins)
}

# Create bins
summary_by_length$bin <- create_bins(summary_by_length$length_2)


calculate_summary <- function(data) {
  # Ensure 'data' includes 'bin', 'phase', and 'mean_score' columns
  
  # Calculate Mean, SD, and N directly, avoiding the complication of collapsing into a list
  aggregated_data <- aggregate(mean_score ~ bin + phase, data, FUN = function(x) c(mean = mean(x, na.rm = TRUE)))
  aggregated_data$sd <- aggregate(mean_score ~ bin + phase, data, FUN = sd, na.rm = TRUE)$mean_score
  aggregated_data$n <- aggregate(mean_score ~ bin + phase, data, FUN = length)$mean_score
  
  # Calculate SEM (Standard Error of the Mean)
  aggregated_data$sem <- aggregated_data$sd / sqrt(aggregated_data$n)
  
  # Calculate the 95% CI for the Mean
  z_score <- 1.96  # Z-score for 95% CI
  aggregated_data$ci_lower <- aggregated_data$mean - z_score * aggregated_data$sem
  aggregated_data$ci_upper <- aggregated_data$mean + z_score * aggregated_data$sem
  
  # Fix column names (since 'aggregate' wraps results in a list)
  aggregated_data$mean <- sapply(aggregated_data$mean, `[`, 1)
  
  return(aggregated_data)
}

# Calculate summary statistics
summary_stats <- calculate_summary(summary_by_length)

# Output the result
print(summary_stats)

### Histogram to explore the distribution of different overlap lengths
# Filter the data for length_2 equal to 6
data_for_hist <- filtered_length_4 %>% 
  filter(length_2 == 25) # change the number here to select different overlap length

# Plot histogram of the filtered data
hist(data_for_hist$normalized_alignment_score_2, breaks = 10)

# Generate a combined histogram using ggplot2
ggplot(filtered_length_4, aes(x = normalized_alignment_score_2)) +
  geom_histogram(binwidth = 0.1, fill = "blue", color = "black", alpha = 0.7) +
  facet_wrap(~length_2, scales = "free_y") + # Facet by length_2 to get separate histograms for each length
  labs(title = "Distribution of Normalized Alignment Scores by Length",
       x = "Normalized Alignment Score",
       y = "Frequency") +
  theme_minimal()


#######################################
# Plotting unique alignment scores by length
######################################


# Filtering the data frame to focus on rows where 'length_2' is greater than 4
count_df <- df

# Now, filtering the data frame to focus on rows where 'length_2" is equal to "length_actual", to
# remove those predictions where the sequence length is not equal, which would impact the alignment score
count_df_filtered <- subset(count_df, length_2 == length_actual)

# Further filtering to exclude rows where 'converted_overlap' is 'none'
count_df_filtered <- subset(count_df_filtered, converted_overlap != "none")

# Aggregate to count unique normalized_alignment_score_2 for each length_2
aggregated_data <- count_df_filtered %>%
  group_by(length_2) %>%
  summarise(count = n_distinct(normalized_alignment_score_2))

# Assuming 'aggregated_data' is already created from previous steps
aggregated_data$phase <- (6 - aggregated_data$length_2) %% 3

# Create the plot with colored bars based on phase
plot <- ggplot(aggregated_data, aes(x = length_2, y = count, fill = as.factor(phase))) +
  geom_col(color = "black") +  # Adding black outline to the bars
  scale_x_continuous(breaks = seq(4, 102, by = 4)) +  # Adjusting x-axis limits
  labs(x = "Overlap Length", y = "Unique Alignment Score Counts",
       fill = "Phase") +
  scale_y_continuous(limits = c(0, 30), breaks = seq(0, 30, by = 2), expand = c(0, 0)) +  # Setting y-axis limits from 0 to 30 with breaks by 2
  theme_minimal() +
  theme(
    legend.position = "top",  # Hiding the legend
    axis.text.x = element_text(angle = 45, hjust = 1.0, size = 18),
    axis.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22, vjust = 2),
    axis.title.y.right = element_text(vjust = 1.7),
    plot.title = element_text(size = 18, hjust = 0.5),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 18),
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Adding black border with transparent fill
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.15, "cm"),  # Negative value to place ticks outside
    axis.ticks.margin = unit(0.5, "cm")
  )

# Display the plot
print(plot)
