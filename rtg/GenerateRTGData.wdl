version 1.0

workflow GenerateRTGData {

  input {
    File truth_fasta
    File baseline_vcf
    Array[File] sample_vcfs
 
    Int? preemptible_tries
  }

  call PrepTruthData {
    input:
      truth_fasta=truth_fasta,
      preemptible_tries=preemptible_tries
  }

  call ProcessVcf as ProcessBaselineVcf {
    input:
      sample_vcf=baseline_vcf
  }

  scatter (vcf in sample_vcfs) {
      call ProcessVcf {
        input:
          sample_vcf = vcf,
          preemptible_tries = preemptible_tries
    }
      call EvalVcf {
        input:
          truth_sdf = PrepTruthData.truth_sdf,
          baseline_vcf_gz = ProcessBaselineVcf.vcf_gz,
          sample_vcf_gz = ProcessVcf.vcf_gz,
          preemptible_tries = preemptible_tries
    }
  }

  call GeneratePlots {
    input:
      EvalVcf.weighted_roc_file
  }
}

task PrepTruthData {
  input {
    File truth_fasta
    Int? preemptible_tries
  }

  String basename = basename(basename(truth_fasta, "sta"), "fa")
  Int disk_size = ceil(size(truth_fasta, "GB")) + 20

  command <<<
    rtg format -o ~{basename}.sdf ~{truth_fasta}
  >>>
  
  output {
    File truth_sdf = "~{basename}.sdf"
  }

  runtime {
      docker: "us.gcr.io/broad-dsde-methods/rtg:3.11"
      memory: "10 GB"
      disks: "local-disk " + disk_size + " HDD"
      preemptible: select_first([preemptible_tries, 5])
      cpu: 1
  }
}

task ProcessVcf {
  input {
    File sample_vcf
     # runtime
    Int? preemptible_tries
  }

  Int disk_size = ceil(size(sample_vcf, "GB")) + 20
  
  command <<<
    bgzip ~{sample_vcf}
    tabix -p vcf ~{sample_vcf}.gz
  >>>
  
  output {
    File vcf_gz = "~{sample_vcf}.gz"
  }

  runtime {
      docker: "us.gcr.io/broad-dsde-methods/rtg:3.11"
      memory: "10 GB"
      disks: "local-disk " + disk_size + " HDD"
      preemptible: select_first([preemptible_tries, 5])
      cpu: 1
  }
}

task EvalVcf {
  input {
    File truth_sdf
    File baseline_vcf_gz
    File sample_vcf_gz
     # runtime
    Int? preemptible_tries
  }

  Int disk_size = ceil(size(sample_vcf_gz, "GB") + size(baseline_vcf_gz)) + 20
  String basename = basename(sample_vcf_gz, ".vcf.gz")
  
  command <<<
    rtg vcfeval -b ~{baseline_vcf_gz} -c ~{sample_vcf_gz} -o ~{basename}_vcfeval -t ~{truth_sdf} --squash-ploidy --vcf-score-field=INFO.TLOD
    
  >>>
  
  output {
    File weighted_roc_file = "~{basename}_vcfeval/weighted_roc.tsv.gz"
  }

  runtime {
      docker: "us.gcr.io/broad-dsde-methods/rtg:3.11"
      memory: "10 GB"
      disks: "local-disk " + disk_size + " HDD"
      preemptible: select_first([preemptible_tries, 5])
      cpu: 1
  }
}


task GeneratePlots {
  input {
    Array[File] weighted_roc_files
     # runtime
    Int? preemptible_tries
  }

  
  command <<<
    rtg rocplot ~{weighted_roc_files} --png=rtgrocplot.png --precision-sensitivity
  >>>

  output {
    File rocplot = "rtgrocplot.png"
  }

  runtime {
      docker: "us.gcr.io/broad-dsde-methods/rtg:3.11"
      memory: "10 GB"
      disks: "local-disk " + 100 + " HDD"
      preemptible: select_first([preemptible_tries, 5])
      cpu: 1
  }
}

