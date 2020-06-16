version 1.0

#import "https://api.firecloud.org/ga4gh/v1/tools/andrea_methods:MicrobialAlignmentPipeline/versions/8/plain-WDL/descriptor" as AlignAndMarkDuplicates
import https://raw.githubusercontent.com/broadinstitute/GATK-for-Microbes/ah_initial_wdl/MicrobialAlignmentPipeline.wdl?token=ACZ6TUS37TR3UZXPCMDDFS266DFFM

# import "MicrobialAlignmentPipeline.wdl" as AlignAndMarkDuplicates

workflow MicrobialGenomePipeline {

  meta {
    description: "Takes in a bam or fastq files, align to ref and shifted reference."
  }

  input {
    File input_bam
    File input_bam_index
    String sample_name

    # Full reference is only requred if starting with a CRAM (BAM doesn't need these files)
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File ref_amb
    File ref_ann
    File ref_bwt
    File ref_pac
    File ref_sa
     
    String? m2_extra_args
    String? m2_filter_extra_args

    #Optional runtime arguments
    Int? preemptible_tries
    File? gatk_override
    String? gatk_docker_override

  }

  parameter_meta {
    input_bam: "Full WGS hg38 bam or cram"
    sample_name: "Name of file in final output vcf"
  }

  call ShiftReference {
    input:
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      preemptible_tries = preemptible_tries,
      gatk_override = "gs://broad-dsp-spec-ops/scratch/andrea/gatk-package-shiftfasta.jar"
  }

  call IndexReference as IndexShiftedRef {
    input:
    ref_fasta = ShiftReference.shifted_ref_fasta,
    preemptible_tries = preemptible_tries
  }

  call RevertSam {
    input:
      input_bam = input_bam,
      preemptible_tries = preemptible_tries
  }

  call AlignAndMarkDuplicates.MicrobialAlignmentPipeline as AlignToRef {
    input:
      input_bam = RevertSam.unmapped_bam,
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_amb = ref_amb,
      ref_ann = ref_ann,
      ref_bwt = ref_bwt,
      ref_pac = ref_pac,
      ref_sa = ref_sa,
      preemptible_tries = preemptible_tries
  }

  call AlignAndMarkDuplicates.MicrobialAlignmentPipeline as AlignToShiftedRef {
    input:
      input_bam = RevertSam.unmapped_bam,
      ref_dict = ShiftReference.shifted_ref_dict,
      ref_fasta = ShiftReference.shifted_ref_fasta,
      ref_fasta_index = ShiftReference.shifted_ref_fasta_index,
      ref_amb = IndexShiftedRef.ref_amb,
      ref_ann = IndexShiftedRef.ref_ann,
      ref_bwt = IndexShiftedRef.ref_bwt,
      ref_pac = IndexShiftedRef.ref_pac,
      ref_sa = IndexShiftedRef.ref_sa,
      preemptible_tries = preemptible_tries
  }

  # call CollectWgsMetrics {
  #   input:
  #     input_bam = AlignToMt.mt_aligned_bam,
  #     input_bam_index = AlignToMt.mt_aligned_bai,
  #     ref_fasta = mt_fasta,
  #     ref_fasta_index = mt_fasta_index,
  #     read_length = max_read_length,
  #     coverage_cap = 100000,
  #     preemptible_tries = preemptible_tries
  # }

  call M2 as CallM2 {
    input:
      input_bam = AlignToRef.aligned_bam,
      input_bai = AlignToRef.aligned_bai,
      ref_fasta = ref_fasta,
      ref_fai = ref_fasta_index,
      ref_dict = ref_dict,
      intervals = ShiftReference.unshifted_intervals,
      gatk_override = gatk_override,
      # Everything is called except the control region.
      # m2_extra_args = select_first([m2_extra_args, ""]) + " -L chrM:576-16024 ",
      preemptible_tries = preemptible_tries
  }

  call M2 as CallShiftedM2 {
    input:
      input_bam = AlignToShiftedRef.aligned_bam,
      input_bai = AlignToShiftedRef.aligned_bai,
      ref_fasta = ShiftReference.shifted_ref_fasta,
      ref_fai = ShiftReference.shifted_ref_fasta_index,
      ref_dict = ShiftReference.shifted_ref_dict,
      intervals = ShiftReference.shifted_intervals,
      gatk_override = gatk_override,
      # Everything is called except the control region.
      # m2_extra_args = select_first([m2_extra_args, ""]) + " -L chrM:8025-9144 ",
      preemptible_tries = preemptible_tries
  }

  call LiftoverAndCombineVcfs {
    input:
      shifted_vcf = CallShiftedM2.raw_vcf,
      vcf = CallM2.raw_vcf,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      shiftback_chain = ShiftReference.shiftback_chain,
      preemptible_tries = preemptible_tries
  }

  call MergeStats {
    input:
      shifted_stats = CallShiftedM2.stats,
      non_shifted_stats = CallM2.stats,
      gatk_override = gatk_override,
      preemptible_tries = preemptible_tries
  }

  # This is a temporary task to handle "joint calling" until Mutect2 can produce a GVCF.
  # This proivdes coverage at each base so low coverage sites can be considered ./. rather than 0/0.
  # call CoverageAtEveryBase {
  #   input:
  #     input_bam_regular_ref = AlignToRef.aligned_bam,
  #     input_bam_regular_ref_index = AlignToRef.aligned_bai,
  #     input_bam_shifted_ref = AlignToShiftedRef.aligned_bam,
  #     input_bam_shifted_ref_index = AlignToShiftedRef.aligned_bai,
  #     shiftback_chain = ShiftReference.shiftback_chain,
  #     control_region_shifted_reference_interval_list = control_region_shifted_reference_interval_list,
  #     non_control_region_interval_list = non_control_region_interval_list,
  #     ref_fasta = ref_fasta,
  #     ref_fasta_index = ref_fasta_index,
  #     ref_dict = ref_dict,
  #     shifted_ref_fasta = shifted_ref_fasta,
  #     shifted_ref_fasta_index = shifted_ref_fasta_index,
  #     shifted_ref_dict = shifted_ref_dict
  # }

  output {
    File final_vcf = LiftoverAndCombineVcfs.final_vcf
    File final_vcf_index = LiftoverAndCombineVcfs.final_vcf_index
    File unmapped_bam = RevertSam.unmapped_bam
    File shifted_ref_dict = ShiftReference.shifted_ref_dict
    File shifted_ref_fasta = ShiftReference.shifted_ref_fasta
    File shifted_ref_fasta_index = ShiftReference.shifted_ref_fasta_index
    File shifted_ref_amb = IndexShiftedRef.ref_amb
    File shifted_ref_ann = IndexShiftedRef.ref_ann
    File shifted_ref_bwt = IndexShiftedRef.ref_bwt
    File shifted_ref_pac = IndexShiftedRef.ref_pac
    File shifted_ref_sa = IndexShiftedRef.ref_sa
  }
}

task ShiftReference {
  input {
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    String basename = basename(ref_fasta, ".fasta")

    # runtime
    Int? preemptible_tries
    File? gatk_override
  }

  Int disk_size = ceil(size(ref_fasta, "GB") * 2.5) + 20

  meta {
    description: "Creates a shifted reference file and shiftback chain file"
  }
  parameter_meta {
    ref_fasta: {
      localization_optional: true
    }
    ref_fasta_index: {
      localization_optional: true
    }    
    ref_dict: {
      localization_optional: true
    }
  }
  command <<<
      set -e

      export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

      gatk --java-options "-Xmx2500m" ShiftFasta \
        -R ~{ref_fasta} \
        -O ~{basename}.shifted.fasta \
        --interval-file-name ~{basename} \
        --shift-back-output ~{basename}.shiftback.chain
  >>>
  runtime {
      docker: "us.gcr.io/broad-gatk/gatk:4.1.7.0"
      memory: "2 GB"
      disks: "local-disk " + disk_size + " HDD"
      preemptible: select_first([preemptible_tries, 5])
      cpu: 2
  }
  output {
    File shifted_ref_fasta = "~{basename}.shifted.fasta"
    File shifted_ref_fasta_index = "~{basename}.shifted.fasta.fai"
    File shifted_ref_dict = "~{basename}.shifted.dict"
    File shiftback_chain = "~{basename}.shiftback.chain"
    File unshifted_intervals = "~{basename}.intervals"
    File shifted_intervals = "~{basename}.shifted.intervals"
  }
}

task IndexReference {
  input {
    File ref_fasta    
    Int? preemptible_tries
  }

  Int disk_size = ceil(size(ref_fasta, "GB") * 2.5) + 20
  String basename = basename(ref_fasta)
  
  command <<<
      set -e
      cp ~{ref_fasta} .
      /usr/gitc/bwa index ~{basename}
      ls -al
      find . -name *.pac -print
  >>>
  runtime {
    preemptible: select_first([preemptible_tries, 5])
    memory: "2 GB"
    disks: "local-disk " + disk_size + " HDD"
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.2-1552931386"
  }

  output {
    File ref_amb = "~{basename}.amb"
    File ref_ann = "~{basename}.ann"
    File ref_bwt = "~{basename}.bwt"
    File ref_pac = "~{basename}.pac"
    File ref_sa = "~{basename}.sa"
  }
}


task RevertSam {
  input {
    File input_bam
    String basename = basename(input_bam, ".bam")

    # runtime
    Int? preemptible_tries
  }
  Int disk_size = ceil(size(input_bam, "GB") * 2.5) + 20

  meta {
    description: "Removes alignment information while retaining recalibrated base qualities and original alignment tags"
  }
  parameter_meta {
    input_bam: "aligned bam"
  }
  command {
    java -Xmx1000m -jar /usr/gitc/picard.jar \
    RevertSam \
    INPUT=~{input_bam} \
    OUTPUT_BY_READGROUP=false \
    OUTPUT=~{basename}.bam \
    VALIDATION_STRINGENCY=LENIENT \
    ATTRIBUTE_TO_CLEAR=FT \
    ATTRIBUTE_TO_CLEAR=CO \
    SORT_ORDER=queryname \
    RESTORE_ORIGINAL_QUALITIES=false
  }
  runtime {
    disks: "local-disk " + disk_size + " HDD"
    memory: "2 GB"
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.2-1552931386"
    preemptible: select_first([preemptible_tries, 5])
  }
  output {
    File unmapped_bam = "~{basename}.bam"
  }
}

# task CoverageAtEveryBase {
#   input {
#     File input_bam_regular_ref
#     File input_bam_regular_ref_index
#     File input_bam_shifted_ref
#     File input_bam_shifted_ref_index
#     File shiftback_chain
#     File control_region_shifted_reference_interval_list
#     File non_control_region_interval_list
#     File ref_fasta
#     File ref_fasta_index
#     File ref_dict
#     File shifted_ref_fasta
#     File shifted_ref_fasta_index
#     File shifted_ref_dict

#     Int? preemptible_tries
#   }
#   Int disk_size = ceil(size(input_bam_regular_ref, "GB") + size(input_bam_shifted_ref, "GB") + size(ref_fasta, "GB") * 2) + 20

#   meta {
#     description: "Remove this hack once there's a GVCF solution."
#   }
#   command <<<
#     set -e

#     java -jar /usr/gitc/picard.jar CollectHsMetrics \
#       I=~{input_bam_regular_ref} \
#       R=~{ref_fasta} \
#       PER_BASE_COVERAGE=non_control_region.tsv \
#       O=non_control_region.metrics \
#       TI=~{non_control_region_interval_list} \
#       BI=~{non_control_region_interval_list} \
#       COVMAX=20000 \
#       SAMPLE_SIZE=1

#     java -jar /usr/gitc/picard.jar CollectHsMetrics \
#       I=~{input_bam_shifted_ref} \
#       R=~{shifted_ref_fasta} \
#       PER_BASE_COVERAGE=control_region_shifted.tsv \
#       O=control_region_shifted.metrics \
#       TI=~{control_region_shifted_reference_interval_list} \
#       BI=~{control_region_shifted_reference_interval_list} \
#       COVMAX=20000 \
#       SAMPLE_SIZE=1

#     R --vanilla <<CODE
#       shift_back = function(x) {
#         if (x < 8570) {
#           return(x + 8000)
#         } else {
#           return (x - 8569)
#         }
#       }

#       control_region_shifted = read.table("control_region_shifted.tsv", header=T)
#       shifted_back = sapply(control_region_shifted[,"pos"], shift_back)
#       control_region_shifted[,"pos"] = shifted_back

#       beginning = subset(control_region_shifted, control_region_shifted[,'pos']<8000)
#       end = subset(control_region_shifted, control_region_shifted[,'pos']>8000)

#       non_control_region = read.table("non_control_region.tsv", header=T)
#       combined_table = rbind(beginning, non_control_region, end)
#       write.table(combined_table, "per_base_coverage.tsv", row.names=F, col.names=T, quote=F, sep="\t")

#     CODE
#   >>>
#   runtime {
#     disks: "local-disk " + disk_size + " HDD"
#     memory: "1200 MB"
#     docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.2-1552931386"
#     preemptible: select_first([preemptible_tries, 5])
#   }
#   output {
#     File table = "per_base_coverage.tsv"
#   }
# }

task LiftoverAndCombineVcfs {
  input {
    File shifted_vcf
    File vcf
    String basename = basename(shifted_vcf, ".vcf")

    File ref_fasta
    File ref_fasta_index
    File ref_dict

    File shiftback_chain

    # runtime
    Int? preemptible_tries
  }

  Float ref_size = size(ref_fasta, "GB") + size(ref_fasta_index, "GB")
  Int disk_size = ceil(size(shifted_vcf, "GB") + ref_size) + 20

  meta {
    description: "Lifts over shifted vcf of control region and combines it with the rest of the chrM calls."
  }
  parameter_meta {
    shifted_vcf: "VCF of control region on shifted reference"
    vcf: "VCF of the rest of chrM on original reference"
    ref_fasta: "Original (not shifted) chrM reference"
    shiftback_chain: "Chain file to lift over from shifted reference to original chrM"
  }
  command<<<
    set -e

    java -jar /usr/gitc/picard.jar LiftoverVcf \
      I=~{shifted_vcf} \
      O=~{basename}.shifted_back.vcf \
      R=~{ref_fasta} \
      CHAIN=~{shiftback_chain} \
      REJECT=~{basename}.rejected.vcf

    java -jar /usr/gitc/picard.jar MergeVcfs \
      I=~{basename}.shifted_back.vcf \
      I=~{vcf} \
      O=~{basename}.final.vcf
    >>>
    runtime {
      disks: "local-disk " + disk_size + " HDD"
      memory: "1 GB"
      docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.2-1552931386"
      preemptible: select_first([preemptible_tries, 5])
    }
    output{
        # rejected_vcf should always be empty
        File rejected_vcf = "~{basename}.rejected.vcf"
        File final_vcf = "~{basename}.final.vcf"
        File final_vcf_index = "~{basename}.final.vcf.idx"
    }
}

task M2 {
  input {
    File ref_fasta
    File ref_fai
    File ref_dict
    File input_bam
    File input_bai
    File intervals
    String? m2_extra_args
    Boolean? make_bamout
    File? gga_vcf
    File? gga_vcf_idx
    File? gatk_override
    # runtime
    Int? mem
    Int? preemptible_tries
  }

  String output_vcf = "raw" + ".vcf"
  String output_vcf_index = output_vcf + ".idx"
  Float ref_size = size(ref_fasta, "GB") + size(ref_fai, "GB")
  Int disk_size = ceil(size(input_bam, "GB") + ref_size) + 20

  # Mem is in units of GB but our command and memory runtime values are in MB
  Int machine_mem = if defined(mem) then mem * 1000 else 3500
  Int command_mem = machine_mem - 500

  meta {
    description: "Mutect2 for calling Snps and Indels"
  }
  parameter_meta {
    input_bam: "Aligned Bam"
    gga_vcf: "VCF for force-calling mode"
  }
  command <<<
      set -e

      export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

      # We need to create these files regardless, even if they stay empty
      # TODO change --mitochondria-mode to --microbial-mode
      touch bamout.bam

      gatk --java-options "-Xmx~{command_mem}m" Mutect2 \
        -R ~{ref_fasta} \
        -I ~{input_bam} \
        ~{"--alleles " + gga_vcf} \
        -O ~{output_vcf} \
        -L ~{intervals} \
        ~{true='--bam-output bamout.bam' false='' make_bamout} \
        ~{m2_extra_args} \
        --annotation StrandBiasBySample \
        --mitochondria-mode \
        --max-reads-per-alignment-start 75 \
        --max-mnp-distance 0
  >>>
  runtime {
      docker: "us.gcr.io/broad-gatk/gatk:4.1.7.0"
      memory: machine_mem + " MB"
      disks: "local-disk " + disk_size + " HDD"
      preemptible: select_first([preemptible_tries, 5])
      cpu: 2
  }
  output {
      File raw_vcf = "~{output_vcf}"
      File raw_vcf_idx = "~{output_vcf_index}"
      File stats = "~{output_vcf}.stats"
      File output_bamOut = "bamout.bam"
  }
}

task MergeStats {
  input {
    File shifted_stats
    File non_shifted_stats
    Int? preemptible_tries
    File? gatk_override
  }

  command{
    set -e

    export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

    gatk MergeMutectStats --stats ~{shifted_stats} --stats ~{non_shifted_stats} -O raw.combined.stats
  }
  output {
    File stats = "raw.combined.stats"
  }
  runtime {
      docker: "us.gcr.io/broad-gatk/gatk:4.1.7.0"
      memory: "3 MB"
      disks: "local-disk 20 HDD"
      preemptible: select_first([preemptible_tries, 5])
  }
}
