version 1.0

import "https://raw.githubusercontent.com/broadinstitute/GATK-for-Microbes/ah_initial_wdl/MicrobialAlignmentPipeline.wdl" as AlignAndMarkDuplicates
import "https://api.firecloud.org/ga4gh/v1/tools/jakec:SamToFastq/versions/8/plain-WDL/descriptor" as SamToFastq
#import "https://raw.githubusercontent.com/broadinstitute/GATK-for-Microbes/ah_initial_wdl/SamToFastq.wdl" as SamToFastq
import "https://raw.githubusercontent.com/gatk-workflows/seq-format-conversion/Paired-FASTQ-to-Unmapped-BAM:3.0.0" as FastqToUnmappedBam

workflow MicrobialGenomePipeline {

  meta {
    description: "Takes in a bam or fastq files, aligns to ref and shifted reference."
  }

  input {
    String sample_name
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File ref_amb
    File ref_ann
    File ref_bwt
    File ref_pac
    File ref_sa
     
    # inputs required when starting from bam file
    File? input_bam
    File? input_bam_index

    # iputs required when starting from fastq files
    File? input_fastq1
    File? input_fastq2
    String? readgroup_name
    String? library_name
    String? platform_unit
    String? run_date
    String? platform_name
    String? sequencing_center
    File? fastqunpaired

    Int? num_dangling_bases
    String? m2_extra_args
    String? m2_filter_extra_args
    Boolean? make_bamout
  

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

# only for bam input
  if (defined(input_bam)) {
    call RevertSam {
      input:
        input_bam = input_bam,
        preemptible_tries = preemptible_tries
    }

    call SamToFastq.SamToFastqTest as SamToFastq {
      input:
        inputBam = input_bam,
        sampleName = sample_name,
        memoryGb = 4,
        diskSpaceGb = 100 # TODO see if we can do computations on the input_bam size here
    }
  }

  if (defined(input_fastq1) && defined(input_fastq2)) {
    call FastqToUnmappedBam.ConvertPairedFastQsToUnmappedBamWf as FastqToUnmappedBam {
      input:
        sample_name = sample_name,
        fastq_1 = input_fastq1,
        fastq_2 = input_fastq1,
        readgroup_name = readgroup_name,
        library_name = library_name,
        platform_unit = platform_unit,
        run_date = run_date,
        platform_name = platform_name,
        sequencing_center = sequencing_center,
        gatk_path = "gatk",
        docker = gatk_docker_override
    }
  }

File? fastq1 = select_first([input_fastq1, SamToFastq.fastq1])
File? fastq2 = select_first([input_fastq1, SamToFastq.fastq2])
File? ubam = select_first([RevertSam.unmapped_bam, FastqToUnmappedBam.output_unmapped_bam])
Int num_dangling_bases = select_first([num_dangling_bases, 3])

# pass in 2 fastq files and unmapped bam
  call AlignAndMarkDuplicates.MicrobialAlignmentPipeline as AlignToRef {
    input:
      unmapped_bam = ubam,
      fastq1 = fastq1,
      fastq2 = fastq2,
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
      unmapped_bam = ubam,
      fastq1 = fastq1,
      fastq2 = fastq2,
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

  call M2 as CallM2 {
    input:
      input_bam = AlignToRef.aligned_bam,
      input_bai = AlignToRef.aligned_bai,
      ref_fasta = ref_fasta,
      ref_fai = ref_fasta_index,
      ref_dict = ref_dict,
      intervals = ShiftReference.unshifted_intervals,
      num_dangling_bases = num_dangling_bases,
      make_bamout = make_bamout,
      gatk_override = gatk_override,
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
      num_dangling_bases = num_dangling_bases,
      make_bamout = make_bamout,
      gatk_override = gatk_override,
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

  call Filter {
    input:
      raw_vcf = LiftoverAndCombineVcfs.merged_vcf,
      raw_vcf_index = LiftoverAndCombineVcfs.merged_vcf_index,
      raw_vcf_stats = MergeStats.stats,
      sample_name = sample_name,
      ref_fasta = ref_fasta,
      ref_fai = ref_fasta_index,
      ref_dict = ref_dict,
      gatk_override = gatk_override,
      gatk_docker_override = gatk_docker_override,
      m2_extra_filtering_args = m2_filter_extra_args,
      # vaf_filter_threshold = 0,  # do we need this value?
      preemptible_tries = preemptible_tries
  }

  output {
    File final_vcf = LiftoverAndCombineVcfs.final_vcf
    File final_vcf_index = LiftoverAndCombineVcfs.final_vcf_index
    File filtered_vcf = Filter.filtered_vcf
    File filtered_vcf_idx = Filter.filtered_vcf_idx
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
    description: "Lifts over shifted vcf of the interval region and combines it with the unshifted vcf."
  }
  parameter_meta {
    shifted_vcf: "VCF of the shifted interval region on shifted reference"
    vcf: "VCF of the unshifted interval region on original reference"
    ref_fasta: "Original (not shifted) reference"
    shiftback_chain: "Chain file to lift over from shifted reference to original reference"
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
    Int num_dangling_bases
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
        --microbial-mode \
        --num-matching-bases-in-dangling-end-to-recover num_dangling_bases \
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

task Filter {
  input {
    File ref_fasta
    File ref_fai
    File ref_dict
    File raw_vcf
    File raw_vcf_index
    File raw_vcf_stats
    Float? vaf_cutoff
    String sample_name

    String? m2_extra_filtering_args
    Int max_alt_allele_count
    Float? vaf_filter_threshold
    Float? f_score_beta

    Boolean run_contamination
    Float? verifyBamID

    File? gatk_override
    String? gatk_docker_override

  # runtime
    Int? preemptible_tries
  }

  String output_vcf = sub(sample_name, "(0x20 | 0x9 | 0xD | 0xA)+", "_") + if compress then ".vcf.gz" else ".vcf"
  Float ref_size = size(ref_fasta, "GB") + size(ref_fai, "GB")
  Int disk_size = ceil(size(raw_vcf, "GB") + ref_size) + 20
  
  meta {
    description: "Mutect2 Filtering for calling Snps and Indels"
  }
  parameter_meta {
      vaf_filter_threshold: "Hard cutoff for minimum allele fraction. All sites with VAF less than this cutoff will be filtered."
      f_score_beta: "F-Score beta balances the filtering strategy between recall and precision. The relative weight of recall to precision."
  }
  command <<<
      set -e

      export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

      # We need to create these files regardless, even if they stay empty
      touch bamout.bam

      gatk --java-options "-Xmx2500m" FilterMutectCalls -V ~{raw_vcf} \
        -R ~{ref_fasta} \
        -O filtered.vcf \
        --stats ~{raw_vcf_stats} \
        ~{m2_extra_filtering_args} \
        # --max-alt-allele-count ~{max_alt_allele_count} \
        --microbial-mode 
        # ~{"--min-allele-fraction " + vaf_filter_threshold} \
        # ~{"--f-score-beta " + f_score_beta} \
        # ~{"--contamination-estimate " + max_contamination}
  >>>
  runtime {
      docker: select_first([gatk_docker_override, "us.gcr.io/broad-gatk/gatk:4.1.8.0"])
      memory: "4 MB"
      disks: "local-disk " + disk_size + " HDD"
      preemptible: select_first([preemptible_tries, 5])
      cpu: 2
  }
  output {
      File filtered_vcf = "~{output_vcf}"
      File filtered_vcf_idx = "~{output_vcf_index}"
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
