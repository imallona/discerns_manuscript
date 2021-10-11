## This file contains the rules to remove exons from the annotation and to create a reduced GTF file.



##################
## Count exon truth

rule exon_truth:
    conda:
        op.join('envs', 'discerns_env.yaml')
    input:
        gtf = config["gtf"],
        sim_iso_res = "simulation/simulated_data/simulated_reads.sim.isoforms.results",
        fastq1 = "simulation/simulated_data/simulated_reads_1.fq",
        fastq2 = "simulation/simulated_data/simulated_reads_2.fq"
    log:
        rlog = op.join('logs', 'count_exons_truth.log')
    output:
        "simulation/analysis/GRCh37.85_all_exon_truth.txt"
    threads: 
        config["cores"]
    script:
        "../scripts/count_exon_truth.R"



##########################
## Reduce exon annotation

rule reduce_GTF:
    conda:
        op.join('envs', 'discerns_env.yaml')
    input:
        gtf = config["gtf"],
        truth = "simulation/analysis/GRCh37.85_all_exon_truth.txt"
    output:
        expand("simulation/reduce_GTF/removed_{removed_exon}_unique.txt", removed_exon = config["reduced_exons"]),
        list(config["reduced_exons"].values()),
        me = config["reduced_gtf"]["me"],
        exon = config["reduced_gtf"]["exon"],
        me_exon = config["reduced_gtf"]["me_exon"]
    script:
        "../scripts/reduce_GTF_expressed_exons.R"


rule removed_exons_truth:
    conda:
        op.join('envs', 'discerns_env.yaml')
    input:
        gtf = lambda wildcards: config["reduced_exons"][wildcards.removed_exon],
        truth = "simulation/analysis/GRCh37.85_all_exon_truth.txt"
    output:
        outfile = "simulation/analysis/removed_exon_truth/{removed_exon}_truth.txt"
    script:
        "../scripts/novel_exon_truth_table.R"

# ## write summary table with information about all removed me/exons
# rule removed_exons_table:
#     input:
#         removed_gtf = lambda wildcards: config["reduced_exons"][wildcards.removed_exon],
#         truth = "simulation/analysis/GRCh37.85_all_exon_truth.txt",
#         reduced_gtf = lambda wildcards: config["reduced_gtf"][wildcards.removed_exon]
#     output:
#         "simulation/analysis/removed_exon_truth/removed_{removed_exon}_summary_table.txt"
#     script:
#         "scripts/write_removed_exons_table.R"


rule classify_removed_exons:
    conda:
        op.join('envs', 'discerns_env.yaml')
    input:
        removed = "simulation/reduced_GTF/removed_{removed_exon}_unique.txt"
    output:
        outfile = "simulation/reduced_GTF/removed_{removed_exon}_unique_classified.txt"
    script:
        "../scripts/classify_removed_exons.R"


rule mapped_junction_count:
    conda:
        op.join('envs', 'discerns_env.yaml')
    input:
        removed = "simulation/reduced_GTF/removed_{removed_exon}_unique_classified.txt",
        bam = "simulation/mapping/STAR/{removed_exon}/{test_dirnames}/pass2_Aligned.out_s.bam"
    output:
        outfile = "simulation/analysis/mapped_junction_count/removed_{removed_exon}_unique_classified_{test_dirnames}_junc_count.txt"
    script:
        "../scripts/mapped_junction_count.R"
