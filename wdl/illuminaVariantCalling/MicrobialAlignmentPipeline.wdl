version 1.0

workflow MicrobialAlignmentPipeline {

  meta {
    description: "Uses BWA to align unmapped bam and marks duplicates."
  }

  input {
    File unmapped_bam
    File fastq1
    File fastq2
    String basename = basename(unmapped_bam, ".bam")
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    File ref_amb
    File ref_ann
    File ref_bwt
    File ref_pac
    File ref_sa

    #Optional runtime arguments
    Int? preemptible_tries
  }

  parameter_meta {
    ref_aligned_bam: "Output is aligned duplicate marked coordinate sorted bam."
  }

  call GetBwaVersion

  call AlignAndMarkDuplicates {
    input:
      unmapped_bam = unmapped_bam,
      fastq1 = fastq1,
      fastq2 = fastq2,
      bwa_version = GetBwaVersion.version,
      output_bam_basename = basename + ".realigned",
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      ref_amb = ref_amb,
      ref_ann = ref_ann,
      ref_bwt = ref_bwt,
      ref_pac = ref_pac,
      ref_sa = ref_sa,
      preemptible_tries = preemptible_tries
  }

  output {
    File aligned_bam = AlignAndMarkDuplicates.output_bam
    File aligned_bai = AlignAndMarkDuplicates.output_bam_index
    File duplicate_metrics = AlignAndMarkDuplicates.duplicate_metrics
    File err = AlignAndMarkDuplicates.bwa_stderr_log
  }
}

task AlignAndMarkDuplicates {
  input {
    File unmapped_bam
    File fastq1
    File fastq2
    String bwa_commandline = "bwa mem -K 100000000 -v 3 -t 2 -Y $bash_ref_fasta"
    String bwa_version
    String output_bam_basename
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File ref_amb
    File ref_ann
    File ref_bwt
    File ref_pac
    File ref_sa

    String? read_name_regex

    Int? preemptible_tries
  }

  String basename = basename(unmapped_bam, ".bam")
  String local_ref_fasta = basename(ref_fasta)
  String metrics_filename = basename + ".metrics"
  Float ref_size = size(ref_fasta, "GB") + size(ref_fasta_index, "GB")
  Int disk_size = ceil(size(unmapped_bam, "GB") * 5 + ref_size) + 20

  meta {
    description: "Aligns with BWA and MergeBamAlignment, then Marks Duplicates. Outputs a coordinate sorted bam."
  }
  parameter_meta {
    unmapped_bam: "Unmapped bam"
    bwa_version: "BWA version to be added to header of aligned bam"
  }
  command <<<
    set -o pipefail
    set -e
    #set -x

    cp ~{ref_fasta} .
    cp ~{ref_fasta_index} .
    cp ~{ref_dict} .
    cp ~{ref_amb} .
    cp ~{ref_ann} .
    cp ~{ref_bwt} .
    cp ~{ref_pac} .
    cp ~{ref_sa} .
  
    # set the bash variable needed for the command-line
    # /usr/gitc/bwa index ~{ref_fasta}
     
    bash_ref_fasta=~{local_ref_fasta}
    # echo "ref_fasta=$bash_ref_fasta"
      /usr/gitc/~{bwa_commandline} ~{fastq1} ~{fastq2} 2> >(tee ~{output_bam_basename}.bwa.stderr.log >&2) > ~{basename}.realigned_bam
     java -Xms3000m -jar /usr/gitc/picard.jar \
      MergeBamAlignment \
      VALIDATION_STRINGENCY=SILENT \
      EXPECTED_ORIENTATIONS=FR \
      ATTRIBUTES_TO_RETAIN=X0 \
      ATTRIBUTES_TO_REMOVE=NM \
      ATTRIBUTES_TO_REMOVE=MD \
      ALIGNED_BAM=~{basename}.realigned_bam \
      UNMAPPED_BAM=~{unmapped_bam} \
      OUTPUT=mba.bam \
      REFERENCE_SEQUENCE=~{local_ref_fasta} \
      PAIRED_RUN=true \
      SORT_ORDER="unsorted" \
      IS_BISULFITE_SEQUENCE=false \
      ALIGNED_READS_ONLY=false \
      CLIP_ADAPTERS=false \
      MAX_RECORDS_IN_RAM=2000000 \
      ADD_MATE_CIGAR=true \
      MAX_INSERTIONS_OR_DELETIONS=-1 \
      PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
      PROGRAM_RECORD_ID="bwamem" \
      PROGRAM_GROUP_VERSION="~{bwa_version}" \
      PROGRAM_GROUP_COMMAND_LINE="~{bwa_commandline}" \
      PROGRAM_GROUP_NAME="bwamem" \
      UNMAPPED_READ_STRATEGY=COPY_TO_TAG \
      ALIGNER_PROPER_PAIR_FLAGS=true \
      UNMAP_CONTAMINANT_READS=true \
      ADD_PG_TAG_TO_READS=false

    java -Xms4000m -jar /usr/gitc/picard.jar \
      MarkDuplicates \
      INPUT=mba.bam \
      OUTPUT=md.bam \
      METRICS_FILE=~{metrics_filename} \
      VALIDATION_STRINGENCY=SILENT \
      ~{"READ_NAME_REGEX=" + read_name_regex} \
      OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
      ASSUME_SORT_ORDER="queryname" \
      CLEAR_DT="false" \
      ADD_PG_TAG_TO_READS=false

    java -Xms4000m -jar /usr/gitc/picard.jar \
      SortSam \
      INPUT=md.bam \
      OUTPUT=~{output_bam_basename}.bam \
      SORT_ORDER="coordinate" \
      CREATE_INDEX=true \
      MAX_RECORDS_IN_RAM=300000
  >>>
  runtime {
    preemptible: select_first([preemptible_tries, 5])
    memory: "6 GB"
    cpu: "2"
    disks: "local-disk " + disk_size + " HDD"
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.2-1552931386"
  }
  output {
    File output_bam = "~{output_bam_basename}.bam"
    File output_bam_index = "~{output_bam_basename}.bai"
    File bwa_stderr_log = "~{output_bam_basename}.bwa.stderr.log"
    File duplicate_metrics = "~{metrics_filename}"
  }
}

task GetBwaVersion {
  meta {
    description: "Gets version of BWA"
  }
  command {
    # not setting set -o pipefail here because /bwa has a rc=1 and we dont want to allow rc=1 to succeed because
    # the sed may also fail with that error and that is something we actually want to fail on.
    /usr/gitc/bwa 2>&1 | \
    grep -e '^Version' | \
    sed 's/Version: //'
  }
  runtime {
    memory: "1 GB"
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.2-1552931386"
  }
  output {
    String version = read_string(stdout())
  }
}
