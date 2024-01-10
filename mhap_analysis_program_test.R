source('/amplseq_required_libraries.R')
source('/amplseq_functions.R')

# Create ArgumentParser object
params <- ArgumentParser()
 
# Define arguments
params$add_argument("-c", "--cigar_paths", help="Path to cigar tables")
params$add_argument("-m", "--metadata_source", help="Metadata source file")
params$add_argument("-l", "--ls_locus_remove", help="File with loci to remove") #Move this to a file
params$add_argument("-s", "--sample_id_pattern", help="Sample ID pattern")
params$add_argument("-pc", "--pop_colors", help="Population colors")
params$add_argument("-pl", "--pop_levels", help="Population levels")
params$add_argument("-pm", "--path_to_markers", help="Path to markers file")
params$add_argument("-gn", "--gene_names", help="Gene names")
params$add_argument("-gi", "--gene_ids", help="Gene IDs")
params$add_argument("-ra", "--reference_alleles", help="Reference alleles file")
params$add_argument("-g", "--gff_file", help="Path to GFF file")
params$add_argument("-f", "--fasta_file", help="Path to FASTA file")
params$add_argument("-v", "--variables", help="Variables")
params$add_argument("-cq", "--collection_quarter", help="Collection quarters")
params$add_argument("-fil", "--filters", help="Filters")
params$add_argument("-pp", "--pop_pairs", help="Population pairs")

cigar_paths <- params$cigar_paths
metadata_source <- params$metadata_source
ls_locus_remove <- params$ls_locus_remove
sample_id_pattern <- params$sample_id_pattern

print(cigar_paths)
print(metadata_source)
print(ls_locus_remove)
print(sample_id_pattern)
print(params$pop_colors)
print(params$pop_levels)

pop_colors <- unlist(strsplit(params$pop_colors, ","))
pop_levels <- adjust_camel_case_name(unlist(strsplit(params$pop_levels, ",")))
pop_levels = factor(pop_levels, levels = pop_levels)
path_to_markers <- params$path_to_markers

if ("pop_pairs" %in% names(params)) {
 pop_pairs <- params$pop_pairs
} else {
 pop_pairs = combn(sort(as.character(pop_levels)), 2, function(x) paste(x, collapse = "-"))
}
gene_names <- params$gene_names
gff_file <- params$gff_file
reference_alleles <- params$reference_alleles
fasta_file <- params$fasta_file
gene_ids <- params$gene_ids
variables <- params$variables
collection_quarter <- params$collection_quarter

#Load the cigar objects and metadata and merge them
cigar_object = read_cigar_tables(paths = cigar_paths, sample_id_pattern = sample_id_pattern)
metadata = read.csv(metadata_source)
cigar_object@metadata = merge(cigar_object@metadata,
                               metadata,
                               by = 'Sample_id',
                               all.x = T)

#Remove certain locus if requested
if (file.exists(ls_locus_remove)) {
  loci <- read.table(ls_locus_remove, header = FALSE, colClasses = "character", stringsAsFactors = FALSE)$V1
  for (locus in loci) {
    cigar_object@cigar_table = cigar_object@cigar_table[,!grepl(locus, colnames(cigar_object@cigar_table))]
  }
}

# Define the geographic unit of analysis and assign them colors
if (length(pop_colors) == length(pop_levels)) {
  pop_colors = pop_colors
} else {
  pop_colors = replicate(n, paste0("#", paste0(sample(0:255, 3, replace = TRUE), collapse = "")))
}
names(pop_levels) = pop_colors
 
# Genotyped samples by Municipality over time
# Plot the sample size (number of genotyped samples) by Municipality over time.
# Quarterly plot
dates = sort(unique(cigar_object@metadata$Quarter_of_Collection))
 
plot_temporal_collection_of_samples = cigar_object@metadata %>%
  filter(!is.na(Subnational_level2)) %>% # Remove samples with no geographic information (Controls)
  summarise(nsamples= n(), .by = c(Subnational_level2, Quarter_of_Collection)) %>%
  ggplot(aes(x = Quarter_of_Collection, y = nsamples, fill = factor(Subnational_level2,
                                                                    levels = pop_levels))) +
  geom_col() +
  theme_bw() +
  scale_fill_manual(values = pop_colors) +
  facet_wrap(.~factor(Subnational_level2,
                      levels = pop_levels),
             strip.position = "top", ncol = 5) +
  labs(title = 'Number of Samples collected over time',
       y = "Collected samples by PCD",
       x = "Quarter") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.position =  "none")

pdf("Results/plot_temporal_collection_of_samples_quarterly.pdf", width = 14)
plot_temporal_collection_of_samples
dev.off()
