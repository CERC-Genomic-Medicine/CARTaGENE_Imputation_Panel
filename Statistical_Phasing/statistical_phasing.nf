#!/usr/bin/env nextflow
/*
* AUTHOR: Mohadese Sayahian Dehkordi, <mohadese.sayahiandehkordi@mail.mcgill.ca>
* VERSION: 1.0
* YEAR: 2023
*/

process prep_SNVs {
    errorStrategy 'retry'
    maxRetries 3
    cache "lenient"
    cpus 1
    memory "8GB"
    time "7h"
    input:
    tuple path(snv_vcf), path(snv_index)
    
    output:
    tuple path("*.ref.vcf.gz"), path("*.ref.log")
    
    publishDir "preprocessed_vcfs/", pattern: "*.vcf.gz", mode: "copy"   
    publishDir "preprocessed_vcfs/", pattern: "*.vcf.gz.tbi", mode: "copy"    

    """
    bcftools view -f PASS $snv_vcf | bcftools annotate -x ^INFO/AF,^INFO/AN,^INFO/AC,^FORMAT/GT | bcftools norm -d all -Oz -o ${snv_vcf.getBaseName()}.filtered.vcf.gz
    bcftools tabix --tbi ${snv_vcf.getBaseName()}.filtered.vcf.gz
    """
}

process prep_SVs{
    errorStrategy 'retry'
    maxRetries 3
    cache "lenient"
    cpus 1
    memory "8GB"
    time "7h"
    input:
    tuple path(sv_vcf), path(sv_index)
    
    output:
    tuple path("*.filtered.vcf.gz"), path("*.filtered.vcf.gz.tbi")
    
    publishDir "preprocessed_SV_vcfs/", pattern: "*.vcf.gz", mode: "copy"   
    publishDir "preprocessed_SV_vcfs/", pattern: "*.vcf.gz.tbi", mode: "copy"    

    """
    bcftools view -f PASS $sv_vcf | bcftools annotate -x ^INFO/AF,^INFO/AN,^INFO/AC,^FORMAT/GT | bcftools norm -d all -Oz -o ${sv_vcf.getBaseName()}.filtered.vcf.gz
    bcftools tabix --tbi ${sv_vcf.getBaseName()}.filtered.vcf.gz
    """
}

process get_chr_name_SNVs {
    errorStrategy 'retry'
    maxRetries 3
    cache "lenient"
    cpus 1
    memory "4GB"
    time "1h"
    input:
    tuple path(snv_vcf), path(snv_index)

    output:
    tuple stdout, path(snv_vcf), path(snv_index)

    """
    tabix -l ${snv_vcf} 
    """
}

process get_chr_name_SVs {
    errorStrategy 'retry'
    maxRetries 3
    cache "lenient"
    cpus 1
    memory "4GB"
    time "1h"
    input:
    tuple path(sv_vcf), path(sv_index)

    output:
    tuple stdout, path(sv_vcf), path(sv_index)

    """
    tabix -l ${sv_vcf} 
    """

}

process concat_vcfs {
    errorStrategy 'retry'
    maxRetries 3
    cache "lenient"
    cpus 1
    memory "4GB"
    time "4h"
    input:
    tuple val(chr_name), path(snv_vcf), path(snv_index), path(sv_vcf), path(sv_index)

    output:
    tuple path("*.combined.vcf.gz"), path("*.combined.vcf.gz.tbi")

    publishDir "combined_vcfs/", pattern: "*.vcf.gz", mode: "copy"   
    publishDir "combined_vcfs/", pattern: "*.vcf.gz.tbi", mode: "copy"    

    """
    bcftools concat $snv_vcf $sv_vcf -Oz -o ${snv_vcf.getBaseName()}.combined.vcf.gz
    bcftools index --tbi ${snv_vcf.getBaseName()}.combined.vcf.gz
    """
}
process beagle_statistical_phasing {
    errorStrategy 'retry'
    maxRetries 3
    cache "lenient"
    cpus 8
    memory "32GB"
    time "10h"
    input:
    tuple path(chr), path(index)
    
    output:
    tuple path("*.ref.vcf.gz"), path("*.ref.log")
    
    publishDir "phased/ref_vcfs/", pattern: "*.vcf.gz", mode: "copy"
    publishDir "phased/ref_logs/", pattern: "*.ref.log", mode: "copy"
    
    """
    java -jar -Djava.io.tmpdir=./temp/ -Xmx32g ${params.beagle} window=25.0 overlap=2.5 nthreads=8 gt=$chr out=${chr.getBaseName()}.ref 
    """
}

process remove_singletons {
    cache "lenient"
    cpus 1
    memory "16GB"
    time "5h"
    input:
    tuple path(chr), path(log)
    
    output:
    tuple path("*.vcf.gz"), path("*.vcf.gz.tbi")
    
    publishDir "phased/ref_withoutsingletons_vcfs/", pattern: "*.vcf.gz", mode: "copy"
    publishDir "phased/ref_withoutsingletons_vcfs/", pattern: "*.vcf.gz.tbi", mode: "copy"

    """
    bcftools view $chr -c 2 -Oz -o  ${chr.getBaseName()}.with_out_singletons.vcf.gz
    bcftools index --tbi ${chr.getBaseName()}.with_out_singletons.vcf.gz
    """
}

workflow {

        snv_ch = Channel.fromPath(params.snv_vcf_path).map{ vcf -> [vcf, vcf + ".tbi" ] }
        sv_ch = Channel.fromPath(params.sv_vcf_path).map{ vcf -> [vcf, vcf + ".tbi" ] }

        prep_snv_ch = prep_SNVs(snv_ch)
        prep_sv_ch = prep_SVs(sv_ch)

        snv_with_chr_name_ch = get_chr_name_SNVs(prep_snv_ch)
        sv_with_chr_name_ch = get_chr_name_SVs(prep_sv_ch)

        stat_phasing_ch = snv_with_chr_name_ch.join(sv_with_chr_name_ch)
        stat_phasing_ch_combine = concat_vcfs(stat_phasing_ch)

        phased_vcfs = beagle_statistical_phasing(stat_phasing_ch_combine)

        remove_singletons(phased_vcfs)
}
