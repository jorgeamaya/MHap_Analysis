library(argparse)
library(rmarkdown)

# Create ArgumentParser object
parser <- ArgumentParser()
 
# Define arguments
parser$add_argument("-c", "--cigar_paths", help="Path to cigar tables")
parser$add_argument("-m", "--metadata_source", help="Metadata source file")
parser$add_argument("-l", "--ls_locus_remove", help="File with loci to remove") #Move this to a file
parser$add_argument("-s", "--sample_id_pattern", help="Sample ID pattern")
parser$add_argument("-pc", "--pop_colors", help="Population colors")
parser$add_argument("-pl", "--pop_levels", help="Population levels")
parser$add_argument("-pm", "--path_to_markers", help="Path to markers file")
parser$add_argument("-gn", "--gene_names", help="Gene names")
parser$add_argument("-gi", "--gene_ids", help="Gene IDs")
parser$add_argument("-ra", "--reference_alleles", help="Reference alleles file")
parser$add_argument("-g", "--gff_file", help="Path to GFF file")
parser$add_argument("-f", "--fasta_file", help="Path to FASTA file")
parser$add_argument("-v", "--variables", help="Variables")
parser$add_argument("-cq", "--collection_quarter", help="Collection quarters")
parser$add_argument("-fil", "--filters", help="Filters")
parser$add_argument("-pp", "--pop_pairs", help="Population pairs")

# Parse the command-line arguments
args <- parser$parse_args()

print(args)
  
# Assign variables based on command-line arguments
render("/mhap_analysis_program.Rmd", params = list(
  cigar_paths = args$cigar_paths,
  metadata_source = args$metadata_source,
  ls_locus_remove = args$ls_locus_remove,
  sample_id_pattern = args$sample_id_pattern,
  pop_colors = args$pop_colors,
  pop_levels = args$pop_levels,
  path_to_markers = args$path_to_markers,
  gene_names = args$gene_names,
  gff_file = args$gff_file,
  reference_alleles = args$reference_alleles,
  fasta_file = args$fasta_file,
  gene_ids = args$gene_ids,
  variables = args$variables,
  collection_quarter = args$collection_quarter
),
  out_dir = "/Results/")
