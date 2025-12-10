
SHELL := /bin/bash

VERSION := $(shell cargo metadata --no-deps --format-version 1 | jq -r '.packages[0].version')




build:
	cargo build --release
	

build-musl:
	cargo build --release --target x86_64-unknown-linux-musl

build-docker:
	bash ./build_linux_dist.sh

build-conda: build-targz
	export version=$(VERSION) && \
	conda build conda_recipe/

build-targz: build-docker build build-musl
	tar -zcvf ver_release/isopedia-$(VERSION).musl.tar.gz -C target/x86_64-unknown-linux-musl/release/ isopedia isopedia-tools -C /ssd1/stix-iso-devspace/isopedia-dev/script/ isopedia-splice-viz.py isopedia-splice-viz-temp.html
	tar -zcvf ver_release/isopedia-$(VERSION).linux.tar.gz -C linux_build/ isopedia isopedia-tools -C /ssd1/stix-iso-devspace/isopedia-dev/script/ isopedia-splice-viz.py isopedia-splice-viz-temp.html

merge2main:
	git checkout main && git merge dev && git tag $(VERSION) && git push && git checkout dev

pack: build build-docker build-conda 

debug:
	export RUST_LOG="debug"
nodebug:
	export RUST_LOG="info"


CMT ?= WIP
git:
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


t2:build
	/usr/bin/time -v target/release/isopedia merge  -i test/HG002.manifest.txt -o test/HG002_idx/

t3:build
	target/release/isopedia index  -i test/HG002_idx/ -m test/HG002.manifest.txt

t4:build
	/usr/bin/time -v target/release/isopedia isoform -i test/HG002_idx/ -g test/gencode.v47.basic.annotation.gtf -o test/test.output.gz --asm

t40:build
	/usr/bin/time -v target/release/isopedia isoform -i test/HG002_idx/ -g test/gencode.v47.basic.annotation.mini.gtf  -o test/test.output.mini.gz --asm

t41:build
	target/release/isopedia isoform -i test/HG002_idx/ -g test/isoseq_transcripts.sorted.filtered_lite.gff -o test/test.output2.gz

t42:build
	target/release/isopedia isoform  -i /ssd1/stix-iso-devspace/stix-isoform-experiment/stage/lrgasp/human_merged_idx \
	 -g /ssd1/stix-iso-devspace/stix-isoform-experiment/data/LRGASP/human_simulation/ground_truth/hs_GENCODE38.basic_annotation.gtf \
	 -o test/test.assembled.output.gz

#grep -w -F 'ENST00000670780.1' /ssd1/stix-iso-devspace/stix-isoform-experiment/data/LRGASP/human_simulation/ground_truth/hs_GENCODE38.basic_annotation.gtf  > test/simulated_R9_missing.gtf


t43:build
	target/release/isopedia isoform --asm  -i /ssd1/stix-iso-devspace/stix-isoform-experiment/stage/lrgasp/human_merged_idx \
	 -g test/simulated_R9_missing.gtf \
	 -o test/test.assembled2.output.gz --info

tlarge:build
	/usr/bin/time -v  target/release/isopedia isoform --asm -n 8  -i /ssd2/isopedia_wk/isopedia_index \
	 -g /ssd1/stix-iso-devspace/stix-isoform-experiment/data/LRGASP/human_simulation/ground_truth/hs_GENCODE38.basic_annotation.gtf \
	 -o test/test.em.output2.gz 

tlarge1:build
	/usr/bin/time -v  target/release/isopedia isoform --asm -n 8  -i /ssd2/isopedia_wk/isopedia_index -c 10 \
	 -g test/chr1.gtf \
	 -o test/test.em.output.chr1.gz 


t45:build
	/usr/bin/time -v  target/release/isopedia isoform --asm  -i test/HG002_idx/ \
	 -g test/gencode49_and_chess313_classcodeU.sorted.chrM.GTF  \
	 -o test/test.em.output2.gz 


tbench:build
	/usr/bin/time -v  target/release/isopedia isoform -n 16  -i /ssd1/stix-iso-devspace/stix-isoform-experiment/stage/lrgasp/human_merged_idx \
	 -g /ssd1/stix-iso-devspace/stix-isoform-experiment/data/LRGASP/human_simulation/ground_truth/hs_GENCODE38.basic_annotation.gtf \
	 -o test/test.em.output.gz 


tbench13:build
	/usr/bin/time -v  target/release/isopedia isoform -n 16  -i /ssd1/stix-iso-devspace/stix-isoform-experiment/stage/lrgasp/human_merged_idx \
	 -g test/bench_chr13.gtf  \
	 -o test/test.em.chr13.output.gz 

tchrm:build
	/usr/bin/time -v  target/release/isopedia isoform --asm  -i /ssd2/isopedia_wk/isopedia_index  -c 10\
	 -g test/gencode49_and_chess313_classcodeU.sorted.chrM.GTF \
	 -o test/test.em.output2.gz


tsamd:build
	/usr/bin/time -v  target/release/isopedia isoform --asm  -i /ssd1/stix-iso-devspace/stix-isoform-experiment/stage/lrgasp/human_merged_idx \
	 -g test/lrgasp_sim_ont.SAMD11.gtf  \
	 -o test/test.em.output3.gz --verbose



tfusion:build
	target/release/isopedia fusion -i test/HG002_idx/ -p chr1:181130,chr1:201853853 -f 200 -o test/fusion.output.gz

tfusion2:build
	target/release/isopedia-anno-fusion -i ../stix-isoform-experiment/stage/fusion_idx/ -f 200 -P test/fusion.bed  -o test/fusion.output.gz


tfusion3:build
	target/release/isopedia fusion -i test/HG002_idx/ -f 200 -G test/gencode.v47.basic.annotation.gtf  -o test/fusion_discovery.output.gz

tfusion4:build
	target/release/isopedia-anno-fusion --verbose -i /ssd1/stix-iso-devspace/stix-isoform-experiment/stage/fusion_idx/ \
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


	target/release/isopedia splice -i test/HG002_idx/  -s 17:7675236,chr17:7675993  -o test/splice.output.gz
	python script/isopedia-splice-viz.py  -i test/splice.output.gz -g /ssd1/stix-iso-devspace/stix-isoform-experiment/ref/gencode.v47.basic.annotation.gtf  -t script/temp.html  -o test/isopedia-splice-view
	zcat test/splice.output.splice.gz  >aa
	zcat test/splice.output.gz  >bb

isoquant:
	isoquant.py  --reference ../stix-isoform-experiment/ref/GRCh38.p14.allChr.fa \
	--genedb ../stix-isoform-experiment/ref/gencode.v47.basic.annotation.gtf \
	--bam  test/bams/hg002.pbcluster.bam \
	--data_type pacbio_ccs \
	-o test/isoquant_output/ \


t1g:build

	target/release/isopedia profile \
	-g /ssd1/stix-iso-devspace/isopedia-dev/test/gencode.v47.basic.annotation.gtf  \
	-o test/gencode_gtf_v47_basic.isoform.gz --gid --tid

	target/release/isopedia profile \
	-g /ssd1/stix-iso-devspace/isopedia-dev/test/gencode.v49.annotation.gtf  \
	-o test/gencode_gtf_v49.isoform.gz --gid --tid

# 	target/release/isopedia profile \
# 	-g /ssd1/stix-iso-devspace/isopedia-dev/ENST00000381578.6.gtf \
# 	-o test/ENST00000381578.isoform.gz --gid --tid --rname


# 	target/release/isopedia profile \
# 	-g /ssd1/stix-iso-devspace/isopedia-dev/ENST00000634662.1.gtf \
# 	-o test/ENST00000634662.isoform.gz --gid --tid --rname


t2g:build

	target/release/isopedia merge -i test/gencode_index_manifest2.tsv -o test/gencode_index/

t3g:build

	target/release/isopedia index -i test/gencode_index/ -m test/gencode_index_manifest2.tsv

t4g:build

	target/release/isopedia isoform -i test/gencode_index/  -g test/gencode.v47.basic.annotation.gtf -o test/gencode_gtf_v47.anno.isoform.gz --info


gtf_index:build t1g t2g t3g t4g

missing:
	cargo build --release && target/release/isopedia isoform -i test/gencode_index --info -f 0 -g test/gencode_missingv47.gtf  -o test/test.output.gz


watchmem:
	@BIN=isopedia; \
	echo "Waiting for process '$$BIN' ..."; \
	while true; do \
		PID=$$(pgrep $$BIN); \
		if [ ! -z "$$PID" ]; then \
			break; \
		fi; \
		sleep 1; \
	done; \
	echo "Monitoring PID $$PID ..."; \
	while true; do \
		date '+%Y-%m-%d %H:%M:%S'; \
		ps -p $$PID -o pid,%mem,rss,vsz,comm; \
		echo "------------------------------"; \
		sleep 4; \
	done | tee mem.log


concat:
	cat src/cmd/isoform.rs  src/grouped_tx.rs src/isoformarchive.rs  src/bptree.rs  > aa.rs