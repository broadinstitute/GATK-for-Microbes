version 1.0

import "Structs.wdl"
import "https://raw.githubusercontent.com/broadinstitute/GATK-for-Microbes/bg_adding_canu_wdl/wdl/longReads/CorrectTrimAssemble.wdl" as CTA

workflow AssembleLongReads {

    meta {
    description: "Workflow to run assemble microbial nanopore or pacbio datasets using canu. If both types of longreads provided, then will run an additional combined assembly "
  }
    input {
        String output_file_prefix
        String genome_size
        String longreads_type
        File nanopore_reads
        File pacbio_reads
        Float nanopore_corrected_error_rate
        Float pacbio_corrected_error_rate


    }

    if(longreads_type == "both") {
        call MultiTechCanu {
            input:
                output_file_prefix = "both" + output_file_prefix,
                genome_size = genome_size,
                nanopore_reads = nanopore_reads,
                pacbio_reads = pacbio_reads,

        }
        call CTA.CorrectTrimAssemble as BothPbCTA {
            input:
                output_file_prefix = "pb" + output_file_prefix,
                genome_size = genome_size,
                reads = pacbio_reads,
                corrected_error_rate = pacbio_corrected_error_rate,
                longreads_type = "pacbio"
        }
        
        call CTA.CorrectTrimAssemble as BothNpCTA {
            input:
                output_file_prefix = "np" + output_file_prefix,
                genome_size = genome_size,
                reads = nanopore_reads,
                corrected_error_rate = nanopore_corrected_error_rate,
                longreads_type = "nanopore"
        }
        
    }

    if(longreads_type == "pacbio") {
        call CTA.CorrectTrimAssemble as PbCTA {
            input:
                output_file_prefix = "pb" + output_file_prefix,
                genome_size = genome_size,
                reads = pacbio_reads,
                corrected_error_rate = pacbio_corrected_error_rate,
                longreads_type = longreads_type
        }
    }

    if(longreads_type == "nanopore") {
        call CTA.CorrectTrimAssemble as NpCTA {
            input:
                output_file_prefix = "np" + output_file_prefix,
                genome_size = genome_size,
                reads = nanopore_reads,
                corrected_error_rate = nanopore_corrected_error_rate,
                longreads_type = longreads_type
        }
    }



    output {
        File? both_contigs_fasta = MultiTechCanu.assembled_contigs_fasta
        File? both_pacbio_contigs_fasta = BothPbCTA.canu_contigs_fasta
        File? both_nanopore_contigs_fasta = BothNpCTA.canu_contigs_fasta
        File? pacbio_contigs_fasta = PbCTA.canu_contigs_fasta
        File? nanopore_contigs_fasta = NpCTA.canu_contigs_fasta
    }
}

# performs assembly on both nanopre and pacbio raw reads
task MultiTechCanu {
    input {
        String output_file_prefix
        String genome_size
        File nanopore_reads
        File pacbio_reads

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        pacbio_reads:           "input raw pacbio reads"
        nanopore_reads:         "input raw nanopore reads"
        genome_size:            "genome size you’re assembling"
        output_file_prefix:     "prefix to output files"
        corrected_error_rate:   "parameter to canu's \'correctedErrorRate\'"
        longreads_type:         "whether the reads are a) pacbio or b) nanopore or c) both"
    }

    Int disk_size = 100 * ceil(size(nanopore_reads, "GB"))

    command <<<
        set -euxo pipefail

        canu -correct \
            -p ~{output_file_prefix} -d both_canu_output \
            genomeSize=~{genome_size} \
            -nanopore ~{nanopore_reads} \
            -pacbio ~{pacbio_reads} \
            || cat /cromwell_root/monitoring.log
    >>>

    output {
        File assembled_contigs_fasta = "both_canu_output/~{output_file_prefix}.contigs.fasta"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          32,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-canu:0.1.0"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}
