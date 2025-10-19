
SHELL := /bin/bash

VERSION := $(shell cargo metadata --no-deps --format-version 1 | jq -r '.packages[0].version')


build:
	cargo build --release
	cargo build --release --target x86_64-unknown-linux-musl

build-docker:
	bash ./build_linux_dist.sh

build-conda: build-targz
	export version=$(VERSION) && \
	conda build conda_recipe/

build-targz: build-docker build
	tar -zcvf ver_release/isopedia-$(VERSION).musl.tar.gz -C target/x86_64-unknown-linux-musl/release/ isopedia isopedia-tools -C /ssd1/stix-iso-devspace/isopedia-dev/script/ isopedia-splice-viz.py temp.html
	tar -zcvf ver_release/isopedia-$(VERSION).linux.tar.gz -C linux_build/ isopedia isopedia-tools -C /ssd1/stix-iso-devspace/isopedia-dev/script/ isopedia-splice-viz.py temp.html

pack: build build-docker build-conda 


CMT ?= WIP
push:
	cargo fmt && git add . && git commit -m "$(CMT)" && git push


t1:build

	target/release/isopedia profile \
	-i /ssd1/stix-iso-devspace/stix-isoform-experiment/data/Kinnex-flrna-DATA-Revio-HG002-1/3-ClusterMap/mapped.bam  \
	-o flnc.isoform.out


t11:build
	target/release/isopedia profile \
	-i test/bams/hg002.directrna.B-Lymphocyte.bam  \
	-o  test/ont.isoform.out.gz --rname

	target/release/isopedia profile \
	-i test/bams/hg002.pbcluster.bam \
	-o test/clustered.mapped.isoform.out.gz --rname

	target/release/isopedia profile \
	-i test/bams/hg002.flnc.minimap2.sorted.bam \
	-o  test/flnc.isoform.out.gz --rname

t1111:build
	target/release/isopedia profile \
	-i test/bams/hg002.flnc.minimap2.sorted.bam \
	-o test/flnc.isoform.out -d 


scp2bcm:
	scp -r -oHostKeyAlgorithms=+ssh-rsa  /ssd1/stix-iso-devspace/isopedia-dev/target/x86_64-unknown-linux-musl/release/isopedia-profile  u249633@hgsc-analysis1:/stornext/snfs170/next-gen/Fritz_Production/zhengxc/isopedia-index/SRA/bin
	scp -r -oHostKeyAlgorithms=+ssh-rsa  /ssd1/stix-iso-devspace/isopedia-dev/target/x86_64-unknown-linux-musl/release/isopedia-tool  u249633@hgsc-analysis1:/stornext/snfs170/next-gen/Fritz_Production/zhengxc/isopedia-index/SRA/bin
	scp -r -oHostKeyAlgorithms=+ssh-rsa  /ssd1/stix-iso-devspace/isopedia-dev/target/x86_64-unknown-linux-musl/release/isopedia-merge  u249633@hgsc-analysis1:/stornext/snfs170/next-gen/Fritz_Production/zhengxc/isopedia-index/SRA/bin
	scp -r -oHostKeyAlgorithms=+ssh-rsa  /ssd1/stix-iso-devspace/isopedia-dev/target/x86_64-unknown-linux-musl/release/isopedia-index  u249633@hgsc-analysis1:/stornext/snfs170/next-gen/Fritz_Production/zhengxc/isopedia-index/SRA/bin


t2:build
	/usr/bin/time -v target/release/isopedia merge  -i test/HG002.manifest.txt -o test/HG002_idx/

t3:build
	target/release/isopedia index  -i test/HG002_idx/ -m test/HG002.manifest.txt

t4:build
	target/release/isopedia isoform -i test/HG002_idx/ -f 10 -g test/gencode.v47.basic.annotation.gtf -o test/test.output.gz

t41:build
	target/release/isopedia isoform -i test/HG002_idx/ -g test/isoseq_transcripts.sorted.filtered_lite.gff -o test/test.output2.gz

tfusion:build
	target/release/isopedia fusion -i test/HG002_idx/ -p chr1:181130,chr1:201853853 -f 200 -o test/fusion.output.gz

tfusion2:build
	target/release/isopedia-anno-fusion -i ../stix-isoform-experiment/stage/fusion_idx/ -f 200 -P test/fusion.bed  -o test/fusion.output.gz


tfusion3:build
	target/release/isopedia fusion -i test/HG002_idx/ -f 200 -G test/gencode.v47.basic.annotation.gtf  -o test/fusion_discovery.output.gz

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


	target/release/isopedia splice -i test/HG002_idx/  -s 17:7675236,17:7675993  -o test/splice.output.gz
	python script/isopedia-splice-viz.py  -i test/splice.output.gz -g /ssd1/stix-iso-devspace/stix-isoform-experiment/ref/gencode.v47.basic.annotation.gtf  -t script/temp.html  -o test/isopedia-splice-view
	zcat test/splice.output.splice.gz  >aa
	zcat test/splice.output.gz  >bb

isoquant:
	isoquant.py  --reference ../stix-isoform-experiment/ref/GRCh38.p14.allChr.fa \
	--genedb ../stix-isoform-experiment/ref/gencode.v47.basic.annotation.gtf \
	--bam  test/bams/hg002.pbcluster.bam \
	--data_type pacbio_ccs \
	-o test/isoquant_output


t1g:build

# 	target/release/isopedia profile \
# 	-g /ssd1/stix-iso-devspace/isopedia-dev/test/gencode.v47.basic.annotation.gtf  \
# 	-o test/gencode_gtf_v47_basic.isoform.gz --gid --tid --rname

# 	target/release/isopedia profile \
# 	-g /ssd1/stix-iso-devspace/isopedia-dev/test/gencode.v49.annotation.gtf  \
# 	-o test/gencode_gtf_v49.isoform.gz --gid --tid --rname

	target/release/isopedia profile \
	-g /ssd1/stix-iso-devspace/isopedia-dev/ENST00000381578.6.gtf \
	-o test/ENST00000381578.isoform.gz --gid --tid --rname


	target/release/isopedia profile \
	-g /ssd1/stix-iso-devspace/isopedia-dev/ENST00000634662.1.gtf \
	-o test/ENST00000634662.isoform.gz --gid --tid --rname


t2g:build

	target/release/isopedia merge -i test/gencode_index_manifest2.tsv -o test/gencode_index/

t3g:build

	target/release/isopedia index -i test/gencode_index/ -m test/gencode_index_manifest2.tsv

t4g:build

# 	target/release/isopedia isoform -i test/gencode_index/ -f 10  -g /ssd1/stix-iso-devspace/isopedia-dev/ENST00000634662.1.gtf -o test/ENST00000634662.anno.isoform.gz --info
# 	target/release/isopedia isoform -i test/gencode_index/ -f 10  -g test/gencode.v49.annotation.gtf -o test/gencode_gtf_v49.anno.isoform.gz --info

	target/release/isopedia isoform -i test/gencode_index/ -f 10  -g test/gencode.v47.basic.annotation.gtf -o test/gencode_gtf_v47.anno.isoform.gz --info


gtf_index:build t1g t2g t3g t4g

missing:
	cargo build --release && target/release/isopedia isoform -i test/gencode_index -f 10 -g test/gencode_missingv47.gtf  -o test/test.output.gz