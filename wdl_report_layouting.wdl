version 1.0

workflow report_layouting {
	input {
		String sample_id_pattern = 'SP'
		String pop_colors
		String pop_levels
	}

	call report_layouting_process {
		input:	
			sample_id_pattern = sample_id_pattern ,
			pop_colors = pop_colors,
			pop_levels = pop_levels,
	}

	output {
		File html_report_f = report_layouting_process.html_report
		File plot_temporal_collection_of_samples_quarterly_f = report_layouting_process.plot_temporal_collection_of_samples_quarterly

#		File Distribution_of_COI_by_sampling_location_f = report_layouting_process.Distribution_of_COI_by_sampling_location
#		File Plot_temporal_collection_of_samples_monthly_f = report_layouting_process.Plot_temporal_collection_of_samples_monthly
#		File Distribution_of_number_heterozygous_loci_per_sample_f = report_layouting_process.Distribution_of_number_heterozygous_loci_per_sample
#		File Plot_locus_amplificatin_rate_f = report_layouting_process.Plot_locus_amplificatin_rate
#		File Plot_temporal_collection_of_samples_yearly_f = report_layouting_process.Plot_temporal_collection_of_samples_yearly
#		File Plot_poly_by_pop_over_time_f = report_layouting_process.Plot_poly_by_pop_over_time
#		File Samples_amplification_rate_f = report_layouting_process.Samples_amplification_rate
	}
}

task report_layouting_process {
	input {
		String sample_id_pattern = 'SP'
		String pop_colors
		String pop_levels
	}

	command <<<
	set -euxo pipefail
	#set -x
	mkdir Results

	gsutil ls ~{path_to_cigars}
	gsutil -m cp -r ~{path_to_cigars}* .
	unzip mhap_metadata.zip

	find . -type f
	Rscript Code/mhap_scripts.R -c Data/cigar_tables/ -m ~{metadata_source} -l ~{ls_locus_remove} -s ~{sample_id_pattern} -pc ~{pop_colors} -pl ~{pop_levels} -pm ~{path_to_markers}
	Rscript render_report.R -c cigar_tables/ -m Gates_Colombia_metadata.csv -l locus_remove.csv -s SP -pc ~{pop_colors} -pl ~{pop_levels} -pm markers.csv
	
	>>>

	output {
		File plot_temporal_collection_of_samples_quarterly = "Results/plot_temporal_collection_of_samples_quarterly.pdf"
		File html_report = "mhap_analysis_program.html"

#		File Plot_locus_amplificatin_rate = "File plot_locus_amplificatin_rate.pdf"
#		File Distribution_of_COI_by_sampling_location = "File distribution_of_COI_by_sampling_location.pdf"
#		File Plot_temporal_collection_of_samples_monthly = "File plot_temporal_collection_of_samples_monthly.pdf"
#		File Distribution_of_number_heterozygous_loci_per_sample = "File distribution_of_number_heterozygous_loci_per_sample.pdf"
#		File Plot_temporal_collection_of_samples_yearly = "File plot_temporal_collection_of_samples_yearly.pdf"
#		File Plot_poly_by_pop_over_time = "File plot_poly_by_pop_over_time.pdf"
#		File Samples_amplification_rate = "File samples_amplification_rate.pdf"
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
