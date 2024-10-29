library(ShortRead)
library(dada2)
library(Biostrings)
library(ggplot2)
library(readxl)
library(dplyr)
library(tidyr)

setwd("/Users/bross/Desktop/AIMS/Analysis/NGS")

######## Plot read counts ################################################################################

# Specify the path where your FASTQ files are located
fastq_path <- "/Users/bross/Desktop/AIMS/Data/NGS"
fastq_files <- list.files(fastq_path, pattern = "\\.fastq\\.gz$", full.names = TRUE)

# Initialize a vector to hold the read counts
read_counts <- numeric(length(fastq_files))

# Count reads in each filtered file
for (i in seq_along(fastq_files)) {
  # Load the ShortReadQ object
  short_read_data <- readFastq(fastq_files[i])
  # Count the reads
  read_counts[i] <- length(short_read_data)
}

# Create a data frame for plotting
results_df <- data.frame(
  Sample = basename(fastq_files),  # Use the file names as sample names
  Read_Count = read_counts
)


############ Sobstitute sample names to indexes ################################################################################

plate_data <- read_excel("index_names.xlsx")

# Step 2: Extract forward and reverse primers
fwd_primers <- plate_data[[1]]
rev_primers <- colnames(plate_data)[-1]  # First row (excluding the header)

# Step 3: Create a data frame for combinations
# Extract the samples from the data (excluding the first row and column)
sample_matrix <- as.matrix(plate_data)
samples <- as.vector(sample_matrix)

# Create a data frame of combinations
combination_df <- expand.grid(Fwd = fwd_primers, Rev = rev_primers)

# Add the corresponding sample names to the combinations
samples <- samples[-(1:4)]
combination_df$Sample <- samples
combination_df <- combination_df %>%
  arrange(Fwd)

combination_df_R1_R2 <- combination_df %>%
  mutate(Sample_R1 = paste(Sample, "R1", sep = "_"),  # Create R1 sample names
         Sample_R2 = paste(Sample, "R2", sep = "_"))

combination_df_f <- combination_df_R1_R2 %>% 
  select(Fwd, Rev, Sample = Sample_R1) %>%  # Select only the R1 samples
  bind_rows(select(combination_df_R1_R2, Fwd, Rev, Sample = Sample_R2)) %>% 
  arrange(Fwd, Rev)

results_df$Sample <- combination_df_f$Sample


# Plot the number of reads for each sequence 

p1 <- ggplot(results_df, aes(x = Sample, y = Read_Count)) +
  geom_bar(stat = "identity") +
  coord_flip() +  # Flipping coordinates for better readability
  labs(title = "Number of Reads for Each Sequence", x = "Sample", y = "Read Count") +
  theme_minimal()

ggsave("read_counts.pdf", plot = p1, device="pdf", width = 10, height = 10)



############ SCF049 / SCF055 samples isolation ################################################################################


target_samples <- c("SCF049-1", "SCF049-2", "SCF049-3", "SCF055-1", "SCF055-4", "SCF055-5")

# Assuming combination_df_final is already defined
# Create a regex pattern to match both R1 and R2 variants
pattern <- paste0(paste(target_samples, collapse = "|"), "_R[12]$")

# Filter the combination_df_final to include both R1 and R2 for target samples
filtered_combinations <- combination_df_f %>%
  filter(grepl(pattern, Sample))

# Initialize a vector to hold the matching FASTQ file names
matching_fastq_files <- c()

# Specify the path where your FASTQ files are located
fastq_path <- "/Users/bross/Desktop/AIMS/Data/NGS"
fastq_files <- list.files(fastq_path, pattern = "\\.fastq\\.gz$", full.names = TRUE)

# Loop over each sample to find matching FASTQ files
for (i in 1:nrow(filtered_combinations)) {
  fwd_index <- filtered_combinations$Fwd[i]
  rev_index <- filtered_combinations$Rev[i]
  
  # Format the search pattern for the FASTQ files
  search_pattern <- paste0("Fwd_", fwd_index, "-Rev_", rev_index)
  
  # Look for files matching the pattern for both R1 and R2
  matched_files <- grep(search_pattern, fastq_files, value = TRUE)
  
  # Store the matched files in the list
  matching_fastq_files <- c(matching_fastq_files, matched_files)
}

###### Primer Removal ################################################################################

primer_F <- "GTGACCTATGAACTCAGGAGTCCCGAGCAAAATAACCAAGC"  # Forward primer sequence
primer_R <- "CTGAGACTTGCACATCGCAGCGTTGTGGATAGCACGAAT"  # Reverse primer sequence
cutadapt_path <- "/Users/bross/Library/Python/3.9/bin/cutadapt"
filtered_path <- "/Users/bross/Desktop/AIMS/Analysis/NGS/FilteredReads"


cutadapt_path <- "/Users/bross/Library/Python/3.9/bin/cutadapt"

for (fastq_file in matching_fastq_files) {
  
  # Determine the output file name based on whether it's R1 or R2
  if (grepl("_R1", fastq_file)) {
    # For R1, apply the forward primer
    output_file <- file.path(filtered_path, sub("\\.fastq\\.gz$", "_trimmed.fastq.gz", basename(fastq_file)))
    cmd <- sprintf("%s -g %s -o %s %s", cutadapt_path, primer_F, output_file, fastq_file)
  } else if (grepl("_R2", fastq_file)) {
    # For R2, apply the reverse primer
    output_file <- file.path(filtered_path, sub("\\.fastq\\.gz$", "_trimmed.fastq.gz", basename(fastq_file)))
    cmd <- sprintf("%s -g %s -o %s %s", cutadapt_path, primer_R, output_file, fastq_file)
  }
  
  # Print the command for verification (optional)
  cat("Running command:", cmd, "\n")
  
  # Run the Cutadapt command
  system(cmd)
}


############### Check reads lenght after trimming ################################################################################

trimmed_files <- list.files(filtered_path, pattern = "\\.fastq.gz$", full.names = TRUE)

read_lengths_df <- data.frame(Sample = character(), ReadLength = integer())

# Extract read lengths for each file
for (file in trimmed_files) {
  # Read the trimmed FASTQ file
  fastq_data <- readFastq(file)
  
  # Extract the lengths of each read
  read_lengths <- width(fastq_data)
  
  # Extract sample name from the file path
  sample_name <- tools::file_path_sans_ext(basename(file))
  
  # Store the results in the dataframe
  read_lengths_df <- rbind(read_lengths_df, 
                           data.frame(Sample = sample_name, 
                                      ReadLength = read_lengths))
}

# Check the first few rows of the read lengths dataframe
head(read_lengths_df)

summary(read_lengths_df$ReadLength)

# Plot the distribution of read lengths for each sample using ggplot2
p2 <- ggplot(read_lengths_df, aes(x = ReadLength, fill = Sample)) +
  geom_histogram(binwidth = 1, alpha = 0.7, position = "identity") +
  theme_minimal() +
  xlab("Read Length") +
  ylab("Frequency") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Read Length Distribution per Sample")


ggsave("after_trimming_reads_lenght.pdf", plot = p2, device="pdf", width = 20, height = 10)




######## Quality Score ################################################################################


list_F <- sort(list.files(filtered_path, pattern = "_R1_trimmed.fastq.gz", full.names = TRUE))
list_R <- sort(list.files(filtered_path, pattern = "_R2_trimmed.fastq.gz", full.names = TRUE))

########### Starting with R1

filtered_data <- readFastq(list_F[1])
print(filtered_data)

quality_scores <- quality(filtered_data)

# Check for missing values in the quality scores
summary(quality_scores)
quality_matrix <- as(quality_scores, "matrix")
dim(quality_matrix)
mean_quality <- rowMeans(quality_matrix, na.rm = TRUE)
median_quality <- apply(quality_matrix, 1, median, na.rm = TRUE)

# Summary of the mean quality scores
summary(mean_quality)
summary(median_quality)

#plot quality score per read
pdf("mean_reads_quality_scores_R1.pdf", width = 8, height = 6)
plot(mean_quality, type = "l", xlab = "Cycle", ylab = "Mean Quality Score", main = "Mean Quality Scores Across Cycles - R1")
dev.off()

#plot quality score per cycle
mean_quality_cycle <- colMeans(quality_matrix, na.rm = TRUE)
summary(mean_quality_cycle)

pdf("mean_cycle_quality_scores_R1.pdf", width = 8, height = 6)
plot(mean_quality_cycle, type = "l", 
     xlab = "Cycle", 
     ylab = "Mean Quality Score", 
     main = "Mean Quality Scores Across Cycles - R1")
dev.off()

#############Same for R2

filtered_data <- readFastq(list_R[1])
print(filtered_data)

quality_scores <- quality(filtered_data)

# Check for missing values in the quality scores
summary(quality_scores)
quality_matrix <- as(quality_scores, "matrix")
dim(quality_matrix)
mean_quality <- rowMeans(quality_matrix, na.rm = TRUE)
median_quality <- apply(quality_matrix, 1, median, na.rm = TRUE)

# Summary of the mean quality scores
summary(mean_quality)
summary(median_quality)

pdf("mean_reads_quality_scores_R2.pdf", width = 8, height = 6)
plot(mean_quality, type = "l", xlab = "Cycle", ylab = "Mean Quality Score", main = "Mean Quality Scores Across Cycles - R2")
dev.off()

#plot quality score per cycle
mean_quality_cycle <- colMeans(quality_matrix, na.rm = TRUE)
summary(mean_quality_cycle)

pdf("mean_cycle_quality_scores_R2.pdf", width = 8, height = 6)
plot(mean_quality_cycle, type = "l", 
     xlab = "Cycle", 
     ylab = "Mean Quality Score", 
     main = "Mean Quality Scores Across Cycles - R2")
dev.off()


####################### Quality filtering and trimming

filt_F <- file.path(filtered_path, "dadafiltered", basename(list_F))
filt_R <- file.path(filtered_path, "dadafiltered", basename(list_R))

# Apply quality filtering and trimming
out <- filterAndTrim(list_F, filt_F, list_R, filt_R,
                     truncLen = c(250, 200), # Adjust according to read lengths and quality profile
                     maxN = 0, maxEE = c(2, 2), truncQ = 2, rm.phix = TRUE,
                     compress = TRUE, multithread = TRUE)
total_reads_in <- sum(out[, "reads.in"])
total_reads_out <- sum(out[, "reads.out"])
(total_reads_out / total_reads_in) * 100

######################## Errore rate learning

# Learn error rates for each direction
errF <- learnErrors(filt_F, multithread = TRUE)
errR <- learnErrors(filt_R, multithread = TRUE)

plotErrors(errF, nominalQ = TRUE)

# Dereplicate sequences
derepF <- derepFastq(filt_F, verbose = TRUE)
derepR <- derepFastq(filt_R, verbose = TRUE)

# Perform sequence inference
dadaF <- dada(derepF, err = errF, multithread = TRUE)
dadaR <- dada(derepR, err = errR, multithread = TRUE)

#########################merge reads

mergers <- mergePairs(dadaF, derepF, dadaR, derepR, verbose = TRUE)
seq_table <- makeSequenceTable(mergers)

###################################remove chimeras
seq_table_nochim <- removeBimeraDenovo(seq_table, method = "consensus", multithread = TRUE, verbose = TRUE)

#################### feature table

feature_table <- as.data.frame(t(seq_table_nochim))
feature_table$ASV_code <- paste0("ASV", seq_len(nrow(feature_table)))

colnames(feature_table) <- c("SCF049-1", "SCF049-2", "SCF049-3", "SCF055-1", "SCF055-4", "SCF055-5", "ASV_code")
feature_table <- as.data.frame(feature_table)

feature_table_long <- feature_table %>%
  pivot_longer(cols = c("SCF049-1", "SCF049-2", "SCF049-3", "SCF055-1", "SCF055-4", "SCF055-5"),  # Replace with your actual sequencing codes
               names_to = "Sample",
               values_to = "Count")

write.csv(feature_table, file = "feature_table_DIVs.csv", row.names = TRUE)

#Calculate relative abundance
feature_table_long <- feature_table_long %>%
  group_by(Sample) %>%
  mutate(Relative_Abundance = Count / sum(Count)) %>%
  ungroup()

feature_table_long <- feature_table_long %>%
  group_by(Sample) %>%
  mutate(ASV_code = factor(ASV_code, 
                           levels = ASV_code[order(-Relative_Abundance)])) %>%
  ungroup()

write.csv(feature_table_long, file = "feature_table_long.csv", row.names = TRUE)

set.seed(58584)
P50 <- createPalette(50,  c("#ff0000", "#00ff00", "#0000ff"))
names(P50)<-NULL

all_DIV <- unique(feature_table_long$ASV_code)
color_mapping <- setNames(P50[1:length(all_DIV)], all_DIV)

p3 <- ggplot(feature_table_long, aes(x = Sample, y = Relative_Abundance, fill = ASV_code)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "ASV Profile per Sample",
       x = "Sample",
       y = "Relative Abundance") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  scale_fill_manual(values= color_mapping)  

ggsave("SCF049-055_psbA_profile.pdf", plot = p3, device="pdf", width = 20, height = 10)

############ total amount of ASVs per sample

total_asv_counts <- feature_table_long %>%
  group_by(Sample) %>%
  summarise(Total_ASV_Counts = sum(Count, na.rm = TRUE))
