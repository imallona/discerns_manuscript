#!/usr/bin/env snakemake -s
#
# discerns_manuscript main workflow file
#
# recommended execution:
# 0. install deps using soft_installs.sh
# 1. start the virtenv with `source ~/virtenvs/discerns/bin/activate`
# 2. snakemake -s Snakefile -j 10 -n -p

import os.path as op

configfile: "config.yaml"

# print(config)

WD = config['WD']
GTF = op.basename(config['gtf_url'])
GENOME = op.basename(config['genome_url'])

include: op.join('rules', 'utils.snmk')
include: op.join("rules", "rsem_simulation.smk")
include: op.join('rules', 'further_soft_installs.snmk')
include: "rules/reduce_GTF.smk"
# include: "rules/mapping_comparison.smk"
# include: "rules/mapping.smk"
# include: "rules/predict_novel_splicing_events.smk"
# include: "rules/quantification.smk"
include: "rules/mapping_real_data.smk"

# print(config['gtf'])
# print(config['RSEMREF'])
print(op.join("simulation", config["SAMPLENAME"] + ".isoforms.results"))

rule all:
    input:
        'Rout/R_packages_install_state.txt',
        # config["RSEMREF"] + ".n2g.idx.fa",
        config["FASTQ2"],
        op.join("simulation", config["SAMPLENAME"] + ".isoforms.results"),
        ## managed till here
        op.join("simulation", 'simulated_data', "simulated_reads_1.fq")
        # "simulation/" + config["SAMPLENAME"] + ".isoforms.results"
        # "simulation/" + config["SAMPLENAME"] + ".isoforms.results",
        # "simulation/" + config["SAMPLENAME"] + ".isoforms.results",
         # expand("simulation/simulated_data/simulated_reads_chr19_22_{nr}.fq", nr = [1,2])
        
        # expand("simulation/analysis/removed_exon_truth/{removed_exon}_truth.txt", removed_exon = config["reduced_exons"])
        # "simulation/" + config["SAMPLENAME"] + ".isoforms.results"
        # "simulation/simulated_data/simulated_reads_1.fq"

        
# TODO chunk 1 start ---------------------        
rule get_srr:
    conda:
        op.join('envs', 'discerns_env.yaml')
    output:
        one_un = temp(op.join(WD, 'real_data', config['SRR'] + '_1.fastq')),
        two_un = temp(op.join(WD, 'real_data', config['SRR'] + '_2.fastq')),
        tmp1 = protected(op.join(WD, 'real_data', config['SRR'] + '_1.fastq.gz')),
        tmp2 = protected(op.join(WD, 'real_data', config['SRR'] + '_2.fastq.gz'))
    params:
        srr = config['SRR'],
        path = op.join(WD, 'real_data')
    threads:
        5
    log:
        op.join(WD, 'real_data', config['SRR'] + '_retrieval.log')
    shell:
        """
        mkdir -p {params.path}
        cd {params.path}
        prefetch {params.srr} &> {log}
        fasterq-dump {params.srr} &>> {log}
        pigz --keep {output.one_un}
        pigz --keep {output.two_un}
        """

rule get_GTF:
    conda:
        op.join('envs', 'discerns_env.yaml')
    output:
        op.join(WD, 'annotation', op.basename(GTF))
    params:
        path = op.join(WD, 'annotation'),
        gtf_url = config['gtf_url']
    log:
        op.join('logs', GTF + '_retrieval.log')
    threads:
        1
    shell:
        """
        mkdir -p {params.path}
        
        wget {params.gtf_url} -O {output} &> {log}
        """

rule get_genome:
    conda:
        op.join('envs', 'discerns_env.yaml')
    output:
        gz = op.join(WD, 'genome', op.basename(GENOME)),
        uncomp = temp(op.join(WD, 'genome', op.splitext(op.basename(GENOME))[0])),        
    params:
        path = op.join(WD, 'genome'),
        genome_url = config['genome_url']
    log:
        op.join('logs', 'genome', GTF + '_retrieval.log')
    threads:
        5
    shell:
        """
        mkdir -p {params.path}
        
        wget {params.genome_url} -O {output.gz} &> {log}

        ## chrfy-it
        pigz -p {threads} --uncompress {output.gz}
        sed -i 's/>/>chr/g' {output.uncomp}

        ## in my hands, pigz segfaults (?)
        # pigz -p {threads} --compress --keep {output.uncomp}
        gzip --keep {output.uncomp}
        """
        
## Install the required R packages into the environment     
rule R_env_install:
    conda:
        op.join('envs', 'discerns_env.yaml')
    input:
        script = "scripts/install_R_packages.R"
    output:
        "Rout/R_packages_install_state.txt"
    log:
        "Rout/install_R_packages.Rout"
    conda:
        # "envs/r_scripts.yaml"
        "envs/discerns_env.yaml"
    threads:
        1
    shell:
        '''R CMD BATCH --no-restore --no-save "--args outtxt='{output}' " {input.script} {log}'''



##################

# RSEM

##################

rule run_RSEM_simulation:
    input:
        expand("simulation/simulated_data/simulated_reads_chr19_22_{nr}.fq", nr = [1,2]),
        expand("simulation/{samplename}.stat/{samplename}.model", samplename = config["SAMPLENAME"])
        


##################

# Reduce GTF annotation

##################

rule run_reduce_GTF:
    input:
        expand("simulation/analysis/removed_exon_truth/{removed_exon}_truth.txt", removed_exon = config["reduced_exons"])

rule count_reduced_GTF:
    input:
        expand( "simulation/analysis/mapped_junction_count/removed_{removed_exon}_unique_classified_{test_dirnames}_junc_count.txt", removed_exon = config["reduced_gtf"], test_dirnames = config["star_param"])



##################

# Quality control

##################
rule run_fastqc:
    conda:
        op.join('envs', 'discerns_env.yaml')
    input:
        fastq1 = "simulation/simulated_data/simulated_reads_chr19_22_1.fq",
        fastq2 = "simulation/simulated_data/simulated_reads_chr19_22_2.fq"
    output:
        "simulation/fastqc/simulated_reads_chr19_22_1_fastqc.html",
        "simulation/fastqc/simulated_reads_chr19_22_2_fastqc.html"
    shell:
        "fastqc -o simlation/fastqc/ {input.fastq1} {input.fastq2}"



##################

# Run read mapping

##################

rule run_star:
    input:
        expand("reference/STAR/chr19_22/{which_reduced_gtf}/Genome", which_reduced_gtf = config["reduced_gtf"]),
        expand("simulation/mapping/STAR/{which_reduced_gtf}/{test_dirnames}/pass2_SJ.out.tab", which_reduced_gtf = config["reduced_gtf"], test_dirnames = config["star_param"]),
        expand("simulation/mapping/STAR/{which_reduced_gtf}/{test_dirnames}/pass2_Aligned.out_s.bam.bai", which_reduced_gtf = config["reduced_gtf"], test_dirnames = config["star_param"])

        # expand("simulation/mapping/STAR/{which_reduced_gtf}/{test_dirnames}/{bam_name}_s.bam.bai", which_reduced_gtf = "me_exon",
        # test_dirnames = "outSJfilterDistToOtherSJmin0_outSJfilterOverhangMin6", bam_name = ["Aligned.out", "pass2_Aligned.out"])

rule run_hisat2:
    input:
        expand("reference/hisat2/chr19_22/{which_reduced_gtf}/{which_reduced_gtf}_GRCh37.85_chr19_22.1.ht2", which_reduced_gtf = config["reduced_gtf"]),
        expand("reference/hisat2/chr19_22/{which_reduced_gtf}/{which_reduced_gtf}_GRCh37.85_chr19_22.1.ht2", which_reduced_gtf = config["reduced_gtf"],),
        expand("simulation/mapping/hisat2/{which_reduced_gtf}/{bam_name}_s.bam.bai", which_reduced_gtf = config["reduced_gtf"], bam_name = "hisat2")

rule run_tophat2:
    input:
        expand("simulation/mapping/tophat/{which_reduced_gtf}/{test_dirnames}/accepted_hits_s.bam.bai", which_reduced_gtf = config["reduced_gtf"], test_dirnames = config["tophat_param"])



##################

# Exon quantification

##################

rule run_featureCounts:
    input:
        expand("simulation/quantification/featureCounts/{which_reduced_gtf}/{test_dirnames}/featureCounts.rds", which_reduced_gtf = config["reduced_gtf"], test_dirnames = config["star_param"])

rule run_EQP:
    input:
        expand("simulation/quantification/EQP/{which_reduced_gtf}/{test_dirnames}/pass2_Aligned.out_s-exon.cnt", which_reduced_gtf = config["reduced_gtf"], test_dirnames = config["star_param"])

rule run_Salmon_derived_counts:
    input:
        expand("simulation/analysis/derived_Salmon_counts/{which_reduced_gtf}/{test_dirnames}/salmon_coverage_count.txt", which_reduced_gtf = config["reduced_gtf"], test_dirnames = config["star_param"] )

rule run_Salmon_derived_counts_stringtie:
    input:
        expand("simulation/analysis/stringtie/derived_Salmon_counts/{which_reduced_gtf}/{stringtie_param}/{test_dirnames}/salmon_coverage_count.txt", which_reduced_gtf = "me_exon", stringtie_param = config["stringtie_param"], test_dirnames = "outSJfilterOverhangMin6")



##################

# Mapping comparison

##################

rule mapping_offset_comparison:
    input:
        expand("simulation/analysis/mapped_offset/{prefix}{suffix}", prefix = ["all_reads", "reads_removed_exons"], suffix = ["_read_offset_table.txt", "_offset_without_sc.txt"])

rule quality_scores:
    input:
        expand("simulation/analysis/mapped_offset/sc_quality_score/{mapper}/{which_reduced_gtf}_{test_dirnames}_quality_scores_per_position.pdf", mapper = "STAR", which_reduced_gtf = config["reduced_gtf"], test_dirnames = config["star_param"]),
        expand("simulation/analysis/mapped_offset/sc_quality_score/{mapper}/{which_reduced_gtf}_quality_scores_per_position.pdf", mapper="hisat2", which_reduced_gtf = config["reduced_gtf"]),
        expand("SRR3192428/analysis/sc_quality_score/{mapper}/{which_reduced_gtf}_quality_scores_per_position.pdf", mapper = "hisat2", which_reduced_gtf=config["reduced_gtf"]),
        expand("SRR3192428/analysis/sc_quality_score/{mapper}/{which_reduced_gtf}_{test_dirnames}_quality_scores_per_position.pdf", mapper = "STAR", which_reduced_gtf = config["reduced_gtf"], test_dirnames = "default")

rule mapped_truth_sj:
    input:
        expand("simulation/mapped_truth/{mapper}/{which_reduced_gtf}/hisat2{suffix}", mapper="hisat2", which_reduced_gtf = config["reduced_gtf"], suffix=["_evaluation_SJ_all.txt", "_evaluation_SJ_overl_removed_exons.txt"]),
        expand("simulation/mapped_truth/{mapper}/{which_reduced_gtf}/{test_dirnames}/pass2_Aligned.out{suffix}",  mapper = "STAR", which_reduced_gtf = config["reduced_gtf"], test_dirnames = config["star_param"],  suffix=["_evaluation_SJ_all.txt", "_evaluation_SJ_overl_removed_exons.txt"])

## to delete all result files to rerun the rule
# find . -name *_evaluation_SJ_overl_removed_exons.txt -delete
# find . -name *_evaluation_SJ_all.txt -delete

rule mapped_truth_sj_comparison:
    input:
        expand("simulation/analysis/mapped_sj_eval/{prefix}{suffix}", prefix = ["all_reads", "reads_removed_exons"], suffix = ["_evaluation_SJ_barplot.pdf", "_evaluation_SJ_barplot_percent.pdf", "_BAM_unique_mapped_barplot.pdf", "_evaluation_SJ_barplot_accuracy.pdf", "_evaluation_SJ_ROC.pdf", "_evaluation_SJ_PR.pdf", "_evaluation_SJ_barplot_F1.pdf"])



#################

# Real data

################

rule map_real_data:
    input:
        expand("SRR3192428/mapping/STAR/{which_reduced_gtf}/{test_dirnames}/{bam_name}_chr19_22_s.bam.bai", which_reduced_gtf = config["reduced_gtf"], test_dirnames = "default", bam_name = "pass2_Aligned.out"),
        expand("SRR3192428/mapping/hisat2/{which_reduced_gtf}/hisat2_s.bam.bai", which_reduced_gtf = config["reduced_gtf"] )



#################

# Novel exon prediction

################

rule predict_exons:
    input:
        expand("simulation/analysis/filtered_SJ/{exon_pred_dir}/{which_reduced_gtf}/novel_exons_{test_dirnames}.txt", exon_pred_dir = "package", which_reduced_gtf = config["reduced_gtf"], test_dirnames = config["star_param"])

rule run_extend_gtf:
    input: 
        expand("simulation/reduced_GTF_with_predicted_exons/{exon_pred_dir}/{which_reduced_gtf}/GRCh37.85_chr19_22_novel_exons_{test_dirnames}.gtf", exon_pred_dir = "package", which_reduced_gtf = config["reduced_gtf"], test_dirnames = config["star_param"])


rule make_PR_curves:
    input:
        expand("simulation/analysis/exon_prediction_performance/PR/{exon_pred_dir}/{which_reduced_gtf}/{test_dirnames}/PR_expression.pdf", exon_pred_dir = "package", which_reduced_gtf = config["reduced_gtf"], test_dirnames = config["star_param"])

rule write_predicted_fasta:
    input:
        expand("simulation/transcriptome/{which_reduced_gtf}/GRCh37.85_chr19_22_novel_exons_{test_dirnames}.fasta", which_reduced_gtf = config["reduced_gtf"], test_dirnames = config["star_param"] )

rule run_gffcompare:
    input:
        expand("simulation/analysis/gffcompare/{exon_pred_dir}/{which_reduced_gtf}/{test_dirnames}/prediction.annotated.gtf", exon_pred_dir = "package", which_reduced_gtf = "me_exon", test_dirnames = "outSJfilterOverhangMin6")



#################

# Stringtie evaluation

################

rule run_stringtie:
    input:
        expand("simulation/analysis/stringtie/predictions/{which_reduced_gtf}/{stringtie_param}/novel_exons_{test_dirnames}_stringtie.txt", which_reduced_gtf = config["reduced_gtf"], stringtie_param = config["stringtie_param"], test_dirnames = "outSJfilterOverhangMin6")
            #test_dirnames = config["star_param"])

rule stringtie_PR_curves:
    input:
        expand("simulation/analysis/stringtie/PR/{which_reduced_gtf}/{stringtie_param}/{test_dirnames}/PR_class_expr.pdf", which_reduced_gtf = config["reduced_gtf"], stringtie_param = config["stringtie_param"], test_dirnames = "outSJfilterOverhangMin6")
        #test_dirnames = config["star_param"])

rule run_gffcompare_stringtie:
    input:
        expand("simulation/analysis/stringtie/gffcompare/{which_reduced_gtf}/{stringtie_param}/{test_dirnames}/stringtie.annotated.gtf", which_reduced_gtf = "me_exon", stringtie_param = "minReadCoverage1_minIsoformAbundance0.05", test_dirnames = "outSJfilterOverhangMin6")
