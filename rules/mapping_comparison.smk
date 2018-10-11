## Rules to run the R scripts that create the comparison ggplots.

# Count number of reads that overlap with each exon, 

## But it is hard to know what the correct way to do this is. Count all reads taht overlap with an exon? Count only the reads that were actually simulated from this exon? But then we might discard reads that are wrongly mapped to the exon...


# what about only looking at the novel exon? See how many reads were mapped to the novel exons (only count the reads that are actually simulated from these exons) 
# This way we ignore all annotated exons, but we assume that the mappers perform similarly and this is not the focus of our study. 
# We are only interested in the novel exons, because we want to know if the different mappers manage to map reads across the splice junctions, maybe dependent on the actual expression of the novel exon. If the exon is higher expresed, we expect that the mapping will be much better than for a lowly expressed exon.

# --> plot the true novel exon count together with the exon count of the different mappers?

import os


rule compute_mapped_truth_hisat:
    input:
        bam = "simulation/mapping/{mapper}/{which_reduced_gtf}/{bam_name}_s.bam",
        gtf = config["gtf"],
        sim_iso_res = "simulation/simulated_data/simulated_reads_chr19_22.sim.isoforms.results",
        removed_gtf = lambda wildcards: config["reduced_exons"][wildcards.which_reduced_gtf]
    params:
        outprefix = "simulation/mapped_truth/{mapper}/{which_reduced_gtf}/{bam_name}_"
    conda:
         "../envs/r_scripts.yaml"
    output:
        "simulation/mapped_truth/{mapper}/{which_reduced_gtf}/{bam_name}_mapped_truth.txt",
        "simulation/mapped_truth/{mapper}/{which_reduced_gtf}/{bam_name}_mapped_truth_removed_exons.txt"

    script:
        "../scripts/mapped_truth.R"


rule compute_mapped_truth_star:
    input:
        bam = "simulation/mapping/{mapper}/{which_reduced_gtf}/{test_dirnames}/{bam_name}_s.bam",
        gtf = config["gtf"],
        sim_iso_res = "simulation/simulated_data/simulated_reads_chr19_22.sim.isoforms.results",
        removed_gtf = lambda wildcards: config["reduced_exons"][wildcards.which_reduced_gtf]
    params:
        outprefix = "simulation/mapped_truth/{mapper}/{which_reduced_gtf}/{test_dirnames}/{bam_name}_"
    conda:
         "../envs/r_scripts.yaml"
    output: 
        "simulation/mapped_truth/{mapper}/{which_reduced_gtf}/{test_dirnames}/{bam_name}_mapped_truth.txt",
        "simulation/mapped_truth/{mapper}/{which_reduced_gtf}/{test_dirnames}/{bam_name}_mapped_truth_removed_exons.txt"
    script:
        "../scripts/mapped_truth.R"
