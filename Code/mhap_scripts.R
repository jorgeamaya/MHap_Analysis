source('Code/amplseq_required_libraries.R')
source('Code/amplseq_functions.R')

# Create ArgumentParser object
#parser <- ArgumentParser()

# Define arguments
#parser$add_argument("-p", "--path_to_meta", help="Path to input meta file listing fastqs (required)")
#parser$add_argument("-l", "--ls_locus_remove", help="Loci to remove")
#parser$add_argument("-c", "--cigar_paths", help="Path to cigar tables")
#parser$add_argument("-m", "--metadata_source", help="Metadata source file")
#parser$add_argument("-s", "--sample_id_pattern", help="Sample ID pattern")
#parser$add_argument("-pc", "--pop_colors", help="Population colors")
# parser$add_argument("-pl", "--pop_levels", help="Population levels")
# parser$add_argument("-pm", "--path_to_markers", help="Path to markers file")
# parser$add_argument("-gn", "--gene_names", help="Gene names")
# parser$add_argument("-gi", "--gene_ids", help="Gene IDs")
# parser$add_argument("-ra", "--reference_alleles", help="Reference alleles file")
# parser$add_argument("-g", "--gff_file", help="Path to GFF file")
# parser$add_argument("-f", "--fasta_file", help="Path to FASTA file")
# parser$add_argument("-v", "--variables", help="Variables")
# parser$add_argument("-cq", "--collection_quarter", help="Collection quarters")
# parser$add_argument("-fil", "--filters", help="Filters")
# parser$add_argument("-pp", "--pop_pairs", help="Population pairs")
# 
# # Parse the command-line arguments
# args <- parser$parse_args()
# 
# # Assign variables based on command-line arguments
# path_to_meta <- args$path_to_meta
# ls_locus_remove <- args$ls_locus_remove
# cigar_paths <- args$cigar_paths
# metadata_source <- args$metadata_source
# sample_id_pattern <- args$sample_id_pattern
# pop_colors <- args$pop_colors
# pop_levels <- args$pop_levels
# pop_levels = factor(pop_levels, levels = pop_levels)
# path_to_markers <- args$path_to_markers
# gene_names <- args$gene_names
# gene_ids <- args$gene_ids
# reference_alleles <- args$reference_alleles
# gff_file <- args$gff_file
# fasta_file <- args$fasta_file
# variables <- args$variables
# collection_quarter <- args$collection_quarter
# filters <- args$filters
# pop_pairs <- args$pop_pairs

#Test variables
ls_locus_remove = c("PF3D7_1302900,1G", "PF3D7_0612900,215A", "PvDHFR")
cigar_paths = "cigar_tables/"
metadata_source = 'Gates_Colombia_metadata.csv'
sample_id_pattern = "SP"
pop_colors = c("firebrick3", "dodgerblue3", "gold3", "darkseagreen3", "lightsalmon2")
pop_levels = c("Quibdó", "Buenaventura", "Guapi", "Tumaco", "Puerto Inírida")
pop_levels = factor(pop_levels, levels = pop_levels)
path_to_markers = "markers.csv"
gene_names <- c('PfDHFR', 'PfMDR1', 'PfDHPS', 'PfKelch13', 'PF3D7_1447900')
gene_ids = c('PF3D7_0417200', 'PF3D7_0523000', 'PF3D7_0810800', 'PF3D7_1343700', 'PF3D7_1447900')
reference_alleles = 'drugR_alleles.csv'
gff_file = "reference/3D7/PlasmoDB-59_Pfalciparum3D7.gff"
fasta_file = "reference/3D7/PlasmoDB-59_Pfalciparum3D7_Genome.fasta"
variables = c('Sample_id', 'Subnational_level2', 'Quarter_of_Collection')
collection_quarter = c('2020-Q4','2021-Q1','2021-Q2','2021-Q3','2021-Q4','2022-Q1','2022-Q2','2022-Q3')
filters = c(paste("Subnational_level2;", paste(pop_levels, collapse = ",")),
            paste("Quarter_of_Collection;", paste(collection_quarter, collapse = ",")))
pop_pairs = c("Buenaventura-Guapi", "Buenaventura-Quibdó", "Buenaventura-Tumaco", "Buenaventura-Puerto Inírida",
              "Guapi-Quibdó", "Guapi-Tumaco", "Guapi-Puerto Inírida", "Quibdó-Tumaco", "Puerto Inírida-Quibdó", "Puerto Inírida-Tumaco")

#Load the cigar objects and metadata and merge them
cigar_object = read_cigar_tables(paths = cigar_paths, sample_id_pattern = sample_id_pattern)
metadata = read.csv(metadata_source)
cigar_object@metadata = merge(cigar_object@metadata,
                              metadata,
                              by = 'Sample_id',
                              all.x = T)

# Remove certain locus if requested
if (length(ls_locus_remove) > 0) {
  for (locus in ls_locus_remove) {
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

pdf("Plots/plot_temporal_collection_of_samples_quarterly.pdf")
plot_temporal_collection_of_samples
dev.off()

# Monthly plot
dates = sort(unique(cigar_object@metadata$Month_of_Collection))

plot_temporal_collection_of_samples = cigar_object@metadata %>%
  filter(!is.na(Subnational_level2)) %>% # Remove samples with no geographic information (Controls)
  summarise(nsamples= n(), .by = c(Subnational_level2, Month_of_Collection)) %>%
  ggplot(aes(x = Month_of_Collection, y = nsamples, fill = factor(Subnational_level2,
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

pdf("Plots/plot_temporal_collection_of_samples_monthly.pdf")
plot_temporal_collection_of_samples
dev.off()

# Yearly plot
dates = sort(unique(cigar_object@metadata$Year_of_Collection))

plot_temporal_collection_of_samples = cigar_object@metadata %>%
  filter(!is.na(Subnational_level2)) %>% # Remove samples with no geographic information (Controls)
  summarise(nsamples= n(), .by = c(Subnational_level2, Year_of_Collection)) %>%
  ggplot(aes(x = Year_of_Collection, y = nsamples, fill = factor(Subnational_level2,
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

pdf("Plots/plot_temporal_collection_of_samples_yearly.pdf")
plot_temporal_collection_of_samples
dev.off()

################################################################################
markers = read.csv(path_to_markers)
ampseq = cigar2ampseq(cigar_object, markers = markers, min_abd = 10, min_ratio = .1, remove_controls = T)

ampseq_filtered = locus_amplification_rate(ampseq, threshold = .65)

# Plot Amplification rate per loci ----

plot_locus_amplificatin_rate = ggdraw()+
  draw_plot(ampseq_filtered@plots$amplification_rate_per_locus+
              theme(axis.text = element_text(size = 12),
                    axis.title = element_text(size = 12),
                    legend.text = element_text(size = 12)),
            x = 0,
            y = 0,
            width = 1,
            height = .5)+
  draw_plot(ampseq_filtered@plots$all_loci_amplification_rate+
              theme(axis.text = element_text(size = 12),
                    axis.title = element_text(size = 12)),
            x = 0,
            y = .5,
            width = 1,
            height = .5)

pdf("Plots/plot_locus_amplificatin_rate.pdf")
plot_locus_amplificatin_rate
dev.off()

################################################################################
ampseq_filtered = sample_amplification_rate(ampseq_filtered)

pdf("Plots/samples_amplification_rate.pdf")
ampseq_filtered@plots$samples_amplification_rate+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))
dev.off()

################################################################################
ampseq_filtered = get_polygenomic(ampseq_object = ampseq_filtered, strata = "Subnational_level2", na.rm = FALSE,
                                  filters = NULL)

plot_log_NPolyLoci_by_pop = log_scale_histogram(data = ampseq_filtered@metadata[ampseq_filtered@metadata$NPolyLoci != 0,],
                                                var = "NPolyLoci", binwidth = 1, group_by = "Subnational_level2",
                                                levels = pop_levels, x_label = "Number of heterozygous loci per sample",
                                                fill_color = pop_colors,
                                                y_breaks = c(1, 5,10, 30), ncol = 4)+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12)) +
  labs(title = 'Distribution of number of heterozygous loci per sample', 
       y = "# Samples",
       x = "# heterozygous loci per sample")

pdf("Plots/distribution_of_number_heterozygous_loci_per_sample.pdf")
plot_log_NPolyLoci_by_pop
dev.off()

################################################################################
plot_log_coi_by_pop = log_scale_histogram(data = ampseq_filtered@metadata,
                                          var = "coi",
                                          binwidth = 1,
                                          group_by = "Subnational_level2",
                                          levels = pop_levels, x_label = "COI",
                                          fill_color = pop_colors,
                                          y_breaks = c(1, 10, 100, 400), ncol = 5)+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12)) +
  labs(title = 'Distribution of COI by sampling location',
       y = "# Samples")

pdf("Plots/distribution_of_COI_by_sampling_location.pdf")
plot_log_coi_by_pop
dev.off()

################################################################################
ampseq_filtered@metadata %<>% mutate(Pop_quarter = paste(Subnational_level2, Quarter_of_Collection, sep = '_'))

plot_poly_by_pop_over_time = get_polygenomic(ampseq_object = ampseq_filtered, strata = "Pop_quarter", update_popsummary = F, na.rm = TRUE, filters = NULL)

plot_poly_by_pop_over_time = plot_poly_by_pop_over_time %>% filter(grepl(paste(pop_levels, sep = "|"), pop))%>%
  mutate(
    Population = stringr::str_split(pop, '_', simplify = TRUE)[,1],
    Date = stringr::str_split(pop, '_', simplify = TRUE)[,2],
    prop_poly_lower = case_when(
      prop_poly == 0 ~ 0,
      prop_poly != 0 ~ prop_poly_lower),
    prop_poly_upper = case_when(
      prop_poly == 0 ~ 0,
      prop_poly != 0 ~ prop_poly_upper)
  )%>%
  ggplot(aes(x = Date,
             y = prop_poly,
             ymin = prop_poly_lower,
             ymax = prop_poly_upper,
             fill = factor(Population, levels = pop_levels)))+
  geom_col()+
  geom_errorbar(width = .2)+
  facet_wrap(~factor(Population, levels = pop_levels), ncol = 5)+
  theme_bw()+
  scale_fill_manual(values = pop_colors)+
  labs(title = 'Temporal change of the proportion of polyclonal infections',
       y = "Polyclonal infections",
       x = "Date of Collection")+
  theme(axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 315, vjust = 0),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.position =  "none")

pdf("Plots/plot_poly_by_pop_over_time.pdf")
plot_poly_by_pop_over_time
dev.off()

################################################################################
# generating the loci_object----
loci_object = ampseq2loci(ampseq)

# measuring relatedness----
if(file.exists('col_pairwise_relatedness.csv')){
  pairwise_relatedness = read.csv('col_pairwise_relatedness.csv')
}else{
  pairwise_relatedness = NULL
  
  for(w in 1:500){
    start = Sys.time()
    pairwise_relatedness = rbind(pairwise_relatedness,
                                 pairwise_hmmIBD(loci_object, parallel = T, w = w, n = 500))
    time_diff = Sys.time() - start
    
    print(paste0('step ', w, ' done in ', time_diff, ' secs'))
    
  }
  
  write.csv(pairwise_relatedness, 'col_pairwise_relatedness.csv', row.names = F, quote = F)
}

plot_relatedness_distribution_within = plot_relatedness_distribution(
  pairwise_relatedness = pairwise_relatedness,
  metadata = ampseq@metadata,
  Population = 'Subnational_level2',
  fill_color = pop_colors,
  type_pop_comparison = 'within',
  ncol = 5,
  pop_levels = pop_levels
)

pdf("Plots/plot_relatedness_distribution_within.pdf")
plot_relatedness_distribution_within$plot +
  labs(title = 'Distribution of relatedness by study area')
dev.off()

################################################################################
plot_frac_highly_related_within = plot_frac_highly_related(
  pairwise_relatedness = pairwise_relatedness,
  metadata = ampseq@metadata,
  Population = 'Subnational_level2',
  fill_color = pop_colors,
  threshold = 0.99, type_pop_comparison = 'within', pop_levels = pop_levels)

pdf("Plots/plot_frac_highly_related_within.pdf")
plot_frac_highly_related_within$plot+
  labs(title = 'Highly related')
dev.off()

################################################################################

plot_frac_highly_related_over_quarters_within = plot_frac_highly_related_over_time(pairwise_relatedness = pairwise_relatedness,
                                                                                   metadata = ampseq@metadata,
                                                                                   Population = c('Subnational_level2', 'Month_of_Collection'),
                                                                                   fill_color = pop_colors,
                                                                                   threshold = 0.99,
                                                                                   type_pop_comparison = 'within',
                                                                                   pop_levels = pop_levels, ncol = 1)




pdf("Plots/plot_frac_highly_related_over_quarters_within.pdf")
plot_frac_highly_related_over_quarters_within$plot_frac_highly_related +
  labs(title = 'Proportion of highly related samples over time (IBD > 0.99)',
       y = 'Proportion')
dev.off()

################################################################################

nsamples_vs_highlyR = left_join(plot_frac_highly_related_over_quarters_within$plot_frac_highly_related$data %>% mutate(pop = paste(Pop_comparison, Date_Yi, sep = '_'))%>% ungroup()%>%
                                  select(pop, prop),
                                  metadata %>%
                                  filter(!is.na(Subnational_level2)) %>%
                                  group_by(Subnational_level2, Month_of_Collection) %>%
                                  summarise(nsamples= n()) %>% mutate(pop = paste(Subnational_level2, Month_of_Collection, sep = '_')), by = 'pop')  %>%
                                  filter(prop != 0, Subnational_level2 %in% c(pop_levels))

plot_cor_highlyR_nsamples = ggscatter(nsamples_vs_highlyR, x = "prop", y = "nsamples", 
                                      add = "reg.line", conf.int = TRUE, 
                                      cor.coef = TRUE, cor.method = "spearman",
                                      title= 'Number of collected samples vs. Proportion of highly related samples',
                                      xlab = "Proportion of highly related samples", ylab = "Number of samples",
                                      facet.by = 'Subnational_level2')

pdf("Plots/plot_cor_highlyR_nsamples.pdf")
plot_cor_highlyR_nsamples
dev.off()

################################################################################
pdf("Plots/plot_network.pdf")
plot_network = plot_network(pairwise_relatedness = pairwise_relatedness,
                            threshold = .99,
                            metadata = ampseq_filtered@metadata,
                            sample_id = 'Sample_id',
                            group_by = 'Subnational_level2',
                            levels = pop_levels,
                            colors = pop_colors
)
dev.off()

#REWORK THESE FUNCTIONS ONWARDS
################################################################################

plot_relatedness_distribution_between = plot_relatedness_distribution(
  pairwise_relatedness = pairwise_relatedness,
  metadata = ampseq_filtered@metadata,
  Population = 'Subnational_level2',
  fill_color = rep('gray50', 10),
  type_pop_comparison = 'between',
  ncol = 5,
  pop_levels = c("Buenaventura-Guapi", "Buenaventura-Quibdó", "Buenaventura-Tumaco", "Buenaventura-Puerto Inírida",
                 "Guapi-Quibdó", "Guapi-Tumaco", "Guapi-Puerto Inírida", "Quibdó-Tumaco", "Puerto Inírida-Quibdó", "Puerto Inírida-Tumaco")
)

pdf("Plots/plot_relatedness_distribution_between.pdf")
plot_relatedness_distribution_between$plot+
  labs(title = 'Distribution of relatedness')
dev.off()

################################################################################
plot_frac_highly_related_between = plot_frac_highly_related(
  pairwise_relatedness = pairwise_relatedness,
  metadata = ampseq_filtered@metadata,
  Population = 'Subnational_level2',
  fill_color = rep('gray50', 10),
  threshold = 0.99, type_pop_comparison = 'between', pop_levels = pop_pairs)

pdf("Plots/plot_frac_highly_related_between.pdf")
plot_frac_highly_related_between$plot+
  scale_x_discrete(breaks = pop_pairs,
                   labels = format_pop_pairs(pop_pairs))+
  labs(title = 'Highly related')
dev.off()

################################################################################
ampseq_drug = ampseq_filtered

ampseq_drug@gt = cbind(ampseq_drug@gt,
                       ampseq_drug@discarded_loci$gt[
                         rownames(ampseq_drug@discarded_loci$gt) %in%
                           rownames(ampseq_drug@gt),
                         grepl((concatenate_genes(gene_names)),
                               colnames(ampseq_drug@discarded_loci$gt))]
)

ampseq_drug@markers = rbind(ampseq_drug@markers,
                            ampseq_drug@discarded_loci$markers[
                              grepl((concatenate_genes(gene_names)),
                                    ampseq_drug@discarded_loci$markers$amplicon),])

drug_resistant_haplotypes_plot = drug_resistant_haplotypes(ampseq_drug,
                                                           reference_alleles = reference_alleles,
                                                           gene_names = gene_names,
                                                           gene_ids = gene_ids,
                                                           gff_file = gff_file,
                                                           fasta_file = fasta_file,
                                                           variables = variables, 
                                                           na.var.rm = TRUE,
                                                           filters = filters)

blues = brewer.pal(9, 'Blues')
reds = brewer.pal(9, 'Reds')

pdf("Plots/haplo_freq_plot.pdf", width = 14)
drug_resistant_haplotypes_plot$haplo_freq_plot+
  theme(axis.text.x = element_text(angle = 90))+
  scale_fill_manual(values = c(
    'gray75',blues[6], reds[3], # MDR2
    reds[c(8,7)], blues[3], reds[6], blues[6], reds[c(9,3)], # DFHR
    'gray75',blues[4], reds[5], blues[c(2,6)], reds[c(4,9,6,8)], # DHPS
    blues[6], #K13
    'gray75', blues[c(4,2,3)], reds[c(4,5,7,9)], blues[2], reds[4], blues[c(3,4)], reds[c(5,6,7,9)] # MDR1
  ))+
  labs(y = 'Haplotype prevalence in infected individuals')
dev.off()

################################################################################
# Define the temporal unit of analysis
dates = sort(unique(cigar_object@metadata$Month_of_Collection))

# Define the geographic unite of analysis
#pop_levels = unique(cigar_object@metadata$Subnational_level2)
#pop_levels = pop_levels[!is.na(pop_levels)]
#pop_levels = pop_levels[c(2, 1, 3, 4, 5)]
#pop_colors = c("firebrick3", "dodgerblue3", "gold3", "darkseagreen3", "lightsalmon2")

#plot_temporal_collection_of_samples = cigar_object@metadata %>%
#  filter(!is.na(Subnational_level2)) %>%
#  group_by(Subnational_level2, Month_of_Collection) %>%
#  summarise(nsamples= n())%>%
#  ggplot(aes(x = Month_of_Collection, y = nsamples, fill = factor(Subnational_level2,
#                                                                  levels = pop_levels)))+
#  geom_col()+
#  theme_bw()+
#  scale_fill_manual(values = pop_colors)+
#  facet_wrap(.~factor(Subnational_level2,
#                      levels = pop_levels),
#             strip.position = "top", ncol = 1)+
#  labs(title = 'B) # Samples collected over time',
#       y = "Collected samples by PCD",
#       x = "Months")+
#  theme(axis.text = element_text(size = 12),
#        axis.text.x = element_text(angle = 315, vjust = 0),
#        axis.title = element_text(size = 12),
#        strip.text = element_text(size = 12),
#        legend.position =  "none")+
#  scale_alpha_manual(values = c(1,.6))+
#  scale_x_discrete(limits = dates)

# Get all shape files from LAC countries at national level
# SpatialPolygonsDataFrame at country level (level = 0)
col0.spldf = gadm_sp_loadCountries("COL",level = 0, basefile = "./GADMTools/world.rds.files/rds0/")

# SpatialPolygonsDataFrame at Municipality level (level = 2)
col2.spldf = gadm_sp_loadCountries("COL",level = 2, basefile = "./GADMTools/world.rds.files/rds2/")

# Label for Municipalities
localities = data.frame(Localities = factor(pop_levels,
                                            levels = pop_levels[5:1]),
                        long = -70,
                        lat = c(6, 8, 10, 12, 14),
                        colors = pop_colors[5:1])

# Plot the map

## Polygon dot plot of observed G6PDd frequency by LAC country----

plot_study_sites = ggplot() +
  geom_polygon(data=col0.spldf$spdf, aes(x=long, y=lat, group = group), fill = "floralwhite", color = "gray40", linewidth=.3)

for (i in 1:length(pop_levels)) {
  print(pop_levels[i])
  print(names(pop_levels[i]))
  plot_study_sites = plot_study_sites + 
    geom_polygon(data=col2.spldf[["spdf"]][col2.spldf[["spdf"]][["NAME_2"]] == pop_levels[i],], aes(x=long, y=lat, group = group), fill = names(pop_levels[i]), color = "gray40", linewidth=.3) 
}

#plot_study_sites = plot_study_sites +
#  geom_polygon(data=col2.spldf[["spdf"]][col2.spldf[["spdf"]][["NAME_2"]] == "Quibdó",], aes(x=long, y=lat, group = group), fill = "firebrick3", color = "gray40", linewidth=.3) 
#  geom_polygon(data=col2.spldf[["spdf"]][col2.spldf[["spdf"]][["NAME_2"]] == "Buenaventura",], aes(x=long, y=lat, group = group), fill = "dodgerblue3", color = "gray40", linewidth=.3) +
#  geom_polygon(data=col2.spldf[["spdf"]][col2.spldf[["spdf"]][["NAME_2"]] == "Guapí",], aes(x=long, y=lat, group = group), fill = "gold3", color = "gray40", linewidth=.3) +
#  geom_polygon(data=col2.spldf[["spdf"]][col2.spldf[["spdf"]][["NAME_2"]] == "Tumaco",], aes(x=long, y=lat, group = group), fill = "darkseagreen3", color = "gray40", linewidth=.3) +
#  geom_polygon(data=col2.spldf[["spdf"]][col2.spldf[["spdf"]][["NAME_2"]] == "Puerto Inírida",], aes(x=long, y=lat, group = group), fill = "lightsalmon2", color = "gray40", linewidth=.3) +
#geom_point(data=localities, aes(x=long, y=lat, color=Localities), size = 5) +
#geom_label(data=localities, aes(x=long + c(1.3, 2.3, 1.5, 1.2), y=lat, label = Localities)) +
plot_study_sites = plot_study_sites + 
  theme_bw()+
  labs(title = "Study sites in Colombia",
       x = "Longitude",
       y = "Latitude")+
  scale_x_continuous(limits = c(-80, -65))+
  scale_y_continuous(limits = c(-5, 15))+
  scale_color_manual(values = localities$colors)+
  theme(legend.position = "none",
        text = element_text(size = 10),
        axis.text = element_text(size = 10),
        title = element_text(size = 10))+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))

#Plot study areas
pdf("Plots/plot_study_areas.pdf")
plot_study_sites
dev.off()

################################################################################
plot_poly_by_pop = ampseq_filtered@pop_summary %>% ggplot(aes(x = factor(pop, levels = c(pop_levels, "Total")),
                                                              y = prop_poly,
                                                              fill = factor(pop, levels = c(pop_levels, "Total"))))+
  geom_col(alpha = .6) +
  geom_errorbar(aes(ymin = prop_poly_lower, ymax = prop_poly_upper), width = .2)+
  theme_bw() +
  labs(title = "D) Frequency of polyclonal infections",
       y = "Frecquency") +
  scale_size(range = c(3, 6)) +
  scale_fill_manual(values = c(pop_colors, "gray30"))+
  scale_y_continuous(limits = c(0, 0.3))+
  theme(axis.text = element_text(size = 12),
        axis.title = element_blank(),
        legend.position = "none")

round(plot_poly_by_pop$data %>% filter(pop == 'Total') %>% select(prop_poly) %>% unlist(),3)

################################################################################
nsamples_vs_proppoly = left_join(plot_poly_by_pop_over_time$data %>% select(pop, prop_poly),
                                 metadata %>%
                                   mutate(Quarter_of_Collection = paste(substr(Date_of_Collection, 1, 4), quarters.Date(Date_of_Collection), sep='-'))%>%
                                   group_by(Subnational_level2, Quarter_of_Collection) %>%
                                   summarise(nsamples= n()) %>% mutate(pop = paste(Subnational_level2, Quarter_of_Collection, sep = '_')), by = 'pop') %>%
  filter(Subnational_level2 %in% c(pop_levels))

plot_cor_proppoly_nsamples = ggscatter(nsamples_vs_proppoly, x = "prop_poly", y = "nsamples", 
                                       add = "reg.line", conf.int = TRUE, 
                                       cor.coef = TRUE, cor.method = "spearman",
                                       title = 'F) # Samples vs. Prop. Poly',
                                       xlab = "Proportion of polyclonal samples", ylab = "Number of samples", facet.by = 'Subnational_level2')

pdf("Plots/plot_cor_proppoly_nsamples.pdf")
plot_cor_proppoly_nsamples
dev.off()

################################################################################

plot_log_poly_by_pop_map = ggdraw()+
  draw_plot(plot_poly_by_pop,
            x = .5,
            y = .34,
            width = .5,
            height = .22)+
  draw_plot(plot_log_coi_by_pop,
            x = .5,
            y = .56,
            width = .5,
            height = .22)+
  draw_plot(plot_log_NPolyLoci_by_pop,
            x = .5,
            y = .78,
            width = .5,
            height = .22)+
  draw_plot(plot_poly_by_pop_over_time,
            x = 0,
            y = 0,
            width = .5,
            height = .33)+
  draw_plot(plot_cor_proppoly_nsamples,
            x = .5,
            y = 0,
            width = .5,
            height = .33)

pdf("Plots/plot_log_poly_by_pop_map.pdf")
plot_log_poly_by_pop_map
dev.off()