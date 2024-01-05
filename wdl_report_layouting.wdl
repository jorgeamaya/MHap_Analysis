version 1.0

workflow report_layouting {
	input {
		String path_to_cigars
		String metadata_source
		File ls_locus_remove
		String sample_id_pattern = 'SP'
		String pop_colors
		String pop_levels
		File path_to_markers
	}

	call report_layouting_process {
		input:	
			path_to_cigars = path_to_cigars,
			ls_locus_remove = ls_locus_remove,
			metadata_source = metadata_source,
			sample_id_pattern = sample_id_pattern ,
			pop_colors = pop_colors,
			pop_levels = pop_levels,
			path_to_markers = path_to_markers
	}

	output {
		File Distribution_of_COI_by_sampling_location_f = report_layouting_process.Distribution_of_COI_by_sampling_location
		File Plot_temporal_collection_of_samples_monthly_f = report_layouting_process.Plot_temporal_collection_of_samples_monthly
		File Distribution_of_number_heterozygous_loci_per_sample_f = report_layouting_process.Distribution_of_number_heterozygous_loci_per_sample
		File Plot_temporal_collection_of_samples_quarterly_f = report_layouting_process.Plot_temporal_collection_of_samples_quarterly
		File Plot_locus_amplificatin_rate_f = report_layouting_process.Plot_locus_amplificatin_rate
		File Plot_temporal_collection_of_samples_yearly_f = report_layouting_process.Plot_temporal_collection_of_samples_yearly
		File Plot_poly_by_pop_over_time_f = report_layouting_process.Plot_poly_by_pop_over_time
		File Samples_amplification_rate_f = report_layouting_process.Samples_amplification_rate
	}
}

task report_layouting_process {
	input {
		String path_to_cigars
		String metadata_source
		File ls_locus_remove
		String sample_id_pattern = 'SP'
		String pop_colors
		String pop_levels
		File path_to_markers
	}

	command <<<
	set -euxo pipefail
	#set -x
	mkdir Data

	gsutil ls ~{path_to_cigars}
	gsutil -m cp -r ~{path_to_cigars}* Data/
	find . -type f
	Rscript Code/mhap_scripts.R -c Data/cigar_tables/ -m ~{metadata_source} -l ~{ls_locus_remove} -s ~{sample_id_pattern} -pc ~{pop_colors} -pl ~{pop_levels} -pm ~{path_to_markers}
	>>>

	output {
		File Distribution_of_COI_by_sampling_location = "File distribution_of_COI_by_sampling_location.pdf"
		File Plot_temporal_collection_of_samples_monthly = "File plot_temporal_collection_of_samples_monthly.pdf"
		File Distribution_of_number_heterozygous_loci_per_sample = "File distribution_of_number_heterozygous_loci_per_sample.pdf"
		File Plot_temporal_collection_of_samples_quarterly = "File plot_temporal_collection_of_samples_quarterly.pdf"
		File Plot_locus_amplificatin_rate = "File plot_locus_amplificatin_rate.pdf"
		File Plot_temporal_collection_of_samples_yearly = "File plot_temporal_collection_of_samples_yearly.pdf"
		File Plot_poly_by_pop_over_time = "File plot_poly_by_pop_over_time.pdf"
		File Samples_amplification_rate = "File samples_amplification_rate.pdf"
	}

	runtime {
		cpu: 1
		memory: "15 GiB"
		disks: "local-disk 10 HDD"
		bootDiskSizeGb: 10
		preemptible: 3
		maxRetries: 1
		docker: 'jorgeamaya/report_layouting'
	}
}
