
SHELL := /bin/bash

build:
	cargo build --release
	cargo build --release --target x86_64-unknown-linux-musl

strip:
	cargo strip --target x86_64-unknown-linux-musl
	


push:
	cargo fmt && git add . && git commit -m "WIP" && git push


t1:build

	target/release/isopedia-extr \
	-i /ssd1/stix-iso-devspace/stix-isoform-experiment/data/Kinnex-flrna-DATA-Revio-HG002-1/3-ClusterMap/mapped.bam  \
	-o flnc.isoform.out

t11:build
	target/release/isopedia-extr \
	-i test/bams/hg002.directrna.B-Lymphocyte.bam  \
	-o test/ont.isoform.out

	target/release/isopedia-extr \
	-i test/bams/hg002.pbcluster.bam \
	-o test/clustered.mapped.isoform.out

	target/release/isopedia-extr \
	-i test/bams/hg002.flnc.minimap2.sorted.bam \
	-o test/flnc.isoform.out

t1111:build 
	target/release/isopedia-extr \
	-i test/bams/hg002.flnc.minimap2.sorted.bam \
	-o test/flnc.isoform.out -d 


scp2bcm:
	scp -r -oHostKeyAlgorithms=+ssh-rsa  /ssd1/stix-iso-devspace/isopedia-dev/target/x86_64-unknown-linux-musl/release/isopedia-extr  u249633@hgsc-analysis1:/stornext/snfs170/next-gen/Fritz_Production/zhengxc/isopedia-index/SRA/bin
	scp -r -oHostKeyAlgorithms=+ssh-rsa  /ssd1/stix-iso-devspace/isopedia-dev/target/x86_64-unknown-linux-musl/release/isopedia-tool  u249633@hgsc-analysis1:/stornext/snfs170/next-gen/Fritz_Production/zhengxc/isopedia-index/SRA/bin
	scp -r -oHostKeyAlgorithms=+ssh-rsa  /ssd1/stix-iso-devspace/isopedia-dev/target/x86_64-unknown-linux-musl/release/isopedia-aggr  u249633@hgsc-analysis1:/stornext/snfs170/next-gen/Fritz_Production/zhengxc/isopedia-index/SRA/bin
	scp -r -oHostKeyAlgorithms=+ssh-rsa  /ssd1/stix-iso-devspace/isopedia-dev/target/x86_64-unknown-linux-musl/release/isopedia-idx  u249633@hgsc-analysis1:/stornext/snfs170/next-gen/Fritz_Production/zhengxc/isopedia-index/SRA/bin


t2:build
	/usr/bin/time -v target/release/isopedia-aggr  -i test/HG002.manifest.txt -o test/HG002_idx/

t3:build
	target/release/isopedia-idx  -i test/HG002_idx/ -m test/HG002.manifest.txt

t4:build
	target/release/isopedia-anno-isoform -i test/HG002_idx/ -f 5 -g test/gencode.v47.basic.annotation.gtf -o test/test.output.txt

t41:build
	target/release/isopedia-anno-isoform -i test/HG002_idx/ -f 20 -g test/isoseq_transcripts.sorted.filtered_lite.gff -o test/test.output2.txt

tfusion:build
	target/release/isopedia-anno-fusion -i test/HG002_idx/ -p chr1:181130,chr1:201853853 -f 200 -o test/fusion.output.gz

tfusion2:build
	target/release/isopedia-anno-fusion -i ../stix-isoform-experiment/stage/fusion_idx/ -f 200 -P test/fusion.bed  -o test/fusion.output.gz


tfusion3:build
	target/release/isopedia-anno-fusion -i test/HG002_idx/ -f 200 -G test/gencode.v47.basic.annotation.gtf  -o test/fusion_discovery.output.gz

tfusion4:build
	target/release/isopedia-anno-fusion --debug -i /ssd1/stix-iso-devspace/stix-isoform-experiment/stage/fusion_idx/ \
	-G test/gencode.v47.basic.annotation.RUNX1RUNX1T1.gtf  \
	-o test/fusion_discovery_cancer.output.gz &> aa.log

tfusion5:build
	target/release/isopedia-anno-fusion -i /ssd1/stix-iso-devspace/stix-isoform-experiment/stage/fusion_idx/ \
	-G /ssd1/stix-iso-devspace/stix-isoform-experiment/ref/gencode.v47.basic.annotation.gtf  \
	-o test/fusion_discovery_cancer.output.all.gz &> aa.log



tfusion6:build
	target/release/isopedia-anno-fusion  -i /ssd1/stix-iso-devspace/isopedia-dev/test/HG002_idx \
	-G /ssd1/stix-iso-devspace/stix-isoform-experiment/ref/gencode.v47.basic.annotation.gtf  \
	-o test/hg002_fusion_discovery_cancer.output.all.gz &> aa.log

# chr17   HAVANA  exon    7675994 7676272 .       -       .       gene_id "ENSG00000141510.19"; transcript_id "ENST00000413465.6"; gene_type "protein_coding"; gene_name "TP53"; transcript_t>
# chr17   HAVANA  CDS     7675994 7676272 .       -       0       gene_id "ENSG00000141510.19"; transcript_id "ENST00000413465.6"; gene_type "protein_coding"; gene_name "TP53"; transcript_t>
# chr17   HAVANA  exon    7675053 7675236 .       -       .       gene_id "ENSG00000141510.19"; transcript_id "ENST00000413465.6"; gene_type "protein_coding"; gene_name "TP53"; transcript_t>
# chr17   HAVANA  CDS     7675053 7675236 .       -       0       gene_id "ENSG00000141510.19"; transcript_id "ENST00000413465.6"; gene_type "protein_coding"; gene_name "TP53"; transcript_t>
# chr17   HAVANA  exon    7674859 7674971 .       -       .       gene_id "ENSG00000141510.19"; transcript_id "ENST00000413465.6"; gene_type "protein_coding"; gene_name "TP53"; transcript_t>
# chr17   HAVANA  CDS     7674859 7674971 .       -       2       gene_id "ENSG00000141510.19"; transcript_id "ENST00000413465.6"; gene_type "protein_coding"; gene_name "TP53"; transcript_t>
# chr17   HAVANA  exon    7674181 7674290 .       -       .       gene_id "ENSG00000141510.19"; transcript_id "ENST00000413465.6"; gene_type "protein_coding"; gene_name "TP53"; transcript_t>
# chr17   HAVANA  CDS     7674181 7674290 .       -       0       gene_id "ENSG00000141510.19"; transcript_id "ENST00000413465.6"; gene_type "protein_coding"; gene_name "TP53"; transcript_t>
# chr17   HAVANA  exon    7661779 7662014 .       -       .       gene_id "ENSG00000141510.19"; transcript_id "ENST00000413465.6"; gene_type "protein_coding"; gene_name "TP53"; transcript_t>
# chr17   HAVANA  CDS     7661942 7662014 .       -       1


splice:build


	target/release/isopedia-anno-splice -i  /hdd1/isopedia_datadownload/encode_merge_idx/ -s 17:7675236,17:7675993  -o test/splice.output.gz
	python script/isopedia-splice-viz.py  -i test/splice.output.gz -g /ssd1/stix-iso-devspace/stix-isoform-experiment/ref/gencode.v47.basic.annotation.gtf  -t script/temp.html  -o test/isopedia-splice-view