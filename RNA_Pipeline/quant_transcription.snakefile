#Align reads to the reference genome using STAR
rule star_align:
        input:
                fq1 = "processed/fastq/{sample}.1.fastq",
                fq2 = "processed/fastq/{sample}.2.fastq"
        output:
                bam = "processed/STAR/{sample}/{sample}.Aligned.sortedByCoord.out.bam"
        params:
                prefix = "processed/STAR/{sample}/{sample}.",
                rg = 'ID:1 \"LB:1\tPL:Illumina\tSM:{sample}\tPU:1\"'
        resources:
                mem = 100000
        threads: 8
        shell:
                "/software/CGP/canpipe/live/bin/STAR --runThreadN {threads} --outSAMtype BAM SortedByCoordinate --outWigTyp
e bedGraph "
                "--outWigNorm None --outWigStrand Stranded --outSAMattrRGline {params.rg} "
                "--genomeDir {config[star_index]} --limitBAMsortRAM 76899000000 "
                "--outFileNamePrefix {params.prefix} --readFilesIn {input.fq1} {input.fq2} "

#Index sorted bams
rule index_bams:
        input:
                "processed/STAR/{sample}/{sample}.Aligned.sortedByCoord.out.bam"
        output:
                "processed/STAR/{sample}/{sample}.Aligned.sortedByCoord.out.bam.bai"
        resources:
                mem = 500
        threads: 1
        shell:
                "samtools index {input}"

#Check genotype concordance between RNA-seq and VCF
rule check_genotype_concordance:
        input:
                "processed/STAR/{sample}/{sample}.Aligned.sortedByCoord.out.bam",
                "processed/STAR/{sample}/{sample}.Aligned.sortedByCoord.out.bam.bai"
        output:
                "processed/verifyBamID/{sample}.verifyBamID.bestSM"
        params:
                out_prefix = "processed/verifyBamID/{sample}.verifyBamID"
        resources:
                mem = 2000
        threads: 1
        shell:
                "/nfs/users/nfs_l/lbc/verifyBamID/bin/verifyBamID --vcf {config[vcf_file]} --bam {input} --out {params.out_
prefix} --best --ignoreRG"

#Convert bedgraph to bigwig
rule bedgraph_to_bigwig:
        input:
                "processed/STAR/{sample}/{sample}.Aligned.sortedByCoord.out.bam"
        output:
                bw1 = "processed/bigwig/{sample}.str1.bw",
                bw2 = "processed/bigwig/{sample}.str2.bw"
        params:
                bg1 = "processed/STAR/{sample}/{sample}.Signal.Unique.str1.out.bg",
                bg2 = "processed/STAR/{sample}/{sample}.Signal.Unique.str2.out.bg",
                bg3 = "processed/STAR/{sample}/{sample}.Signal.UniqueMultiple.str1.out.bg",
                bg4 = "processed/STAR/{sample}/{sample}.Signal.UniqueMultiple.str2.out.bg"
        resources:
                mem = 1000
        threads: 1
        shell:
                "{config[bedg2bw_root]} {params.bg1} {config[chromosome_lengths]} {output.bw1} && "
                "{config[bedg2bw_root]} {params.bg2} {config[chromosome_lengths]} {output.bw2} && "
                "rm {params.bg1} && rm {params.bg2} && rm {params.bg3} && rm {params.bg4}"
                
#Sort BAMs by name
rule sort_bam_by_name:
        input:
                "processed/STAR/{sample}/{sample}.Aligned.sortedByCoord.out.bam"
        output:
                temp("processed/sorted_bam/{sample}.Aligned.sortedByName.out.bam")
        threads: 6
        resources:
                mem = 8000
        shell:
                "samtools sort -n -m 1000M -o {output} {input}"

#Quantify expression using featureCounts
rule quantify_featureCounts:
        input:
                bam = "processed/sorted_bam/{sample}.Aligned.sortedByName.out.bam"
        output:
                counts = "processed/featureCounts/{sample}.featureCounts.txt",
                summary = temp("processed/featureCounts/{sample}.featureCounts.txt.summary")
        threads: 1
        resources:
                mem = 1000
        shell:
                "{config[feat_root]} -p -C -D 5000 -d 50 --donotsort -a {config[ensembl_gtf]} -o {output.counts} {input.bam}"

#Make sure that all final output files get created
rule make_all:
        input:
                expand("processed/verifyBamID/{sample}.verifyBamID.bestSM", sample=config["samples"]),
                expand("processed/bigwig/{sample}.str1.bw", sample=config["samples"]),
                expand("processed/featureCounts/{sample}.featureCounts.txt", sample=config["samples"]),
        output:
                "processed/quant_out.txt"
        resources:
                mem = 100
        threads: 1
        shell:
                "echo 'Done' > {output}"


#Make sure that all final output files get created
rule make_macroMap:
        input:
                expand("processed/bigwig/{sample}.str1.bw", sample=config["samples"]),
        output:
                "processed/macroMap/out.txt"
        resources:
                mem = 100
        threads: 1
        shell:
                "echo 'Done' > {output}"
