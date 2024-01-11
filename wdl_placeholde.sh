unzip mhap_metadata.zip

mkdir Results

Rscript /render_report.R -c /mhap_metadata/cigar_tables/ -m /mhap_metadata/Gates_Colombia_metadata.csv -l /mhap_metadata/locus_remove.csv -s SP -pc firebrick3,dodgerblue3,gold3,darkseagreen3,lightsalmon2 -pl Quibdó,Buenaventura,Guapi,Tumaco,PuertoInírid -pm /mhap_metadata/markers.csv
