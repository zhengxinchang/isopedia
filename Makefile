
build:
	cargo build --release


t1:build

	target/release/isopedia-extr \
	-b /ssd1/stix-iso-devspace/stix-isoform-experiment/data/Kinnex-flrna-DATA-Revio-HG002-1/3-ClusterMap/mapped.bam  \
	-o flnc.isoform.out

t11:build
	target/release/isopedia-extr \
	-b test/bams/hg002.directrna.B-Lymphocyte.bam  \
	-o test/ont.isoform.out

	target/release/isopedia-extr \
	-b test/bams/hg002.pbcluster.bam \
	-o test/clustered.mapped.isoform.out

	target/release/isopedia-extr \
	-b test/bams/hg002.flnc.minimap2.sorted.bam \
	-o test/flnc.isoform.out

t1111:build 
	target/release/isopedia-extr \
	-b test/bams/hg002.flnc.minimap2.sorted.bam \
	-o test/flnc.isoform.out -d 



t2:build
	target/release/isopedia-aggr  -i test/HG002.manifest.txt -o test/HG002_idx/

t3:build
	target/release/isopedia-idx  -i test/HG002_idx/


t4:build
	target/release/isopedia-anno-isoform -i test/HG002_idx/ -f 5 -g test/gencode.v47.basic.annotation.gtf -o test/test.output.txt

t41:build
	target/release/isopedia-anno-isoform -i test/HG002_idx/ -f 20 -g test/isoseq_transcripts.sorted.filtered_lite.gtf -o test/test.output2.txt

tfusion:build
	target/release/isopedia-anno-fusion -i test/HG002_idx/ -p chr1:181130,chr1:201853853 -f 200 -o test/fusion.output.txt



t8:build # for flnc reads
	target/release/stix-isoform  extract -b test/bams/hg002.flnc.minimap2.sorted.bam -o test/flnc.isoform.out
	target/release/stix-isoform  extract -b test/bams/hg002.pbcluster.bam -o test/clustered.mapped.isoform.out

t81:build
	target/release/stix-isoform  extract -b test/bams/hg002.directrna.B-Lymphocyte.bam -o test/ont.isoform.out

t9:build
	printf "path\tsample\n" > test/HG002.manifest.txt
	printf "test/clustered.mapped.isoform.out\tclustered_bam\n" >> test/HG002.manifest.txt
	printf "test/flnc.isoform.out\tflnc_bam\n" >> test/HG002.manifest.txt
	printf "test/ont.isoform.out\tont-bam\n" >> test/HG002.manifest.txt
	target/release/stix-isoform  aggr -i test/HG002.manifest.txt -o test/HG002_idx/
	target/release/stix-isoform  idx -i test/HG002_idx/

t91:build
	printf "path\tsample\n" > test/HG002.manifest-ont.txt
	printf "test/ont.isoform.out\tont-bam\n" >> test/HG002.manifest-ont.txt
	target/release/stix-isoform  aggr -i test/HG002.manifest-ont.txt -o test/ont_idx/
	target/release/stix-isoform  idx -i test/ont_idx/

t10:build
	time target/release/stix-isoform  search -i test/HG002_idx/ -g test/isoseq_transcripts.sorted.filtered_lite.gtf -o test/Kinnex.clustered.isoseq_transcripts.report.txt 
	time target/release/stix-isoform  search -i test/HG002_idx/ -g test/gencode.v47.basic.annotation.gtf -o test/Kinnex.clustered.gencode.report.txt
	grep 'transcript' test/isoseq_transcripts.sorted.filtered_lite.gtf | wcl 
	grep 'transcript' test/gencode.v47.basic.annotation.gtf | wcl

t101:build
	time  target/release/stix-isoform  search -i test/Kinnex.clustered_idx/ -g test/isoseq_transcripts.sorted.filtered_lite.gtf -o test/Kinnex.clustered.isoseq_transcripts.report.txt 
	time  target/release/stix-isoform  search -i test/Kinnex.clustered_idx/ -g test/gencode.v47.basic.annotation.gtf -o test/Kinnex.clustered.gencode.report.txt

t102:build
	time  target/release/stix-isoform  search -f 10 -i test/ont_idx/ -g test/isoseq_transcripts.sorted.filtered_lite.gtf -o test/ont.isoseq_transcripts.report.txt


tmp:build

	target/release/stix-isoform  search -i test/HG002_idx/ -f 0  -g test/ENST00000623083.4.gtf -o test/ENST00000623083.4.-f0.report.txt  
	target/release/stix-isoform  search -i test/HG002_idx/ -f 2  -g test/ENST00000623083.4.gtf -o test/ENST00000623083.4.-f2.report.txt  
	target/release/stix-isoform  search -i test/HG002_idx/ -f 4  -g test/ENST00000623083.4.gtf -o test/ENST00000623083.4.-f4.report.txt  
	target/release/stix-isoform  search -i test/HG002_idx/ -f 10  -g test/ENST00000623083.4.gtf -o test/ENST00000623083.4.-f10.report.txt  
	target/release/stix-isoform  search -i test/HG002_idx/ -f 12  -g test/ENST00000623083.4.gtf -o test/ENST00000623083.4.-f12.report.txt  
	grep ".*" test/ENST00000623083.4.-f* |cut -f 1,10,11,12 | grep -v chrom | sed 's/test\/ENST00000623083.4.-f//g' | sed 's/.report.txt:1//' | sort -n