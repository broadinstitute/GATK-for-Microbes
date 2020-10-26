version 1.0

workflow GenerateRTGData {

  input {
    File truth_fasta
    File baseline_vcf
    Array[File] sample_vcfs
 
    Int? preemptible_tries
  }

  call ProcessVcf as ProcessBaselineVcf {
    input:
      sample_vcf=baseline_vcf,
      preemptible_tries = preemptible_tries
  }

  scatter (vcf in sample_vcfs) {
      call ProcessVcf {
        input:
          sample_vcf = vcf,
          preemptible_tries = preemptible_tries
    }
  }

  call EvalVcf {
    input:
      truth_fasta=truth_fasta,
      baseline_vcf_gz = ProcessBaselineVcf.vcf_gz,
      baseline_vcf_gz_index = ProcessBaselineVcf.vcf_gz_index,
      sample_vcf_gz_list = ProcessVcf.vcf_gz,
      sample_vcf_gz_index_list = ProcessVcf.vcf_gz_index,
      preemptible_tries = preemptible_tries
}

  call GeneratePlots {
    input:
      weighted_roc_files = EvalVcf.weighted_roc_files,
      preemptible_tries = preemptible_tries
  }
}


task ProcessVcf {
  input {
    File sample_vcf
     # runtime
    Int? preemptible_tries
  }

  Int disk_size = ceil(size(sample_vcf, "GB")) + 20
  String basename = basename(sample_vcf)
  
  command <<<
    ln -s ~{sample_vcf} ~{basename}
    bgzip ~{basename}
    tabix -p vcf ~{basename}.gz
  >>>
  
  output {
    File vcf_gz = "~{basename}.gz"
    File vcf_gz_index = "~{basename}.gz.tbi"
  }

  runtime {
      docker: "us.gcr.io/broad-dsde-methods/rtg:3.11_091520"
      memory: "5 GB"
      disks: "local-disk " + disk_size + " HDD"
      preemptible: select_first([preemptible_tries, 5])
      cpu: 1
  }
}

task EvalVcf {
  input {
    File truth_fasta
    File baseline_vcf_gz
    File baseline_vcf_gz_index
    Array[File] sample_vcf_gz_list
    Array[File] sample_vcf_gz_index_list
     # runtime
    Int? preemptible_tries
  }

  String truth_basename = basename(basename(truth_fasta, "sta"), ".fa")
  
  command <<<
    rtg format -o ~{truth_basename}.sdf ~{truth_fasta}

    # the seperator below is a literal tab character; \t does not work
    for sample_vcf_gz in ~{sep=' ' sample_vcf_gz_list}  ; do
      # figure out the base name of the file
      sample_basename=`echo $sample_vcf_gz | awk -F '/' '{ print $NF }' | sed -E 's/\.vcf[\.gz]+//'`
      rtg vcfeval -b ~{baseline_vcf_gz} -c $sample_vcf_gz -o ${sample_basename}_vcfeval -t ~{truth_basename}.sdf --squash-ploidy --vcf-score-field=INFO.TLOD
      ln -s ${sample_basename}_vcfeval/weighted_roc.tsv.gz ${sample_basename}_weighted_roc.tsv.gz
    done
    
  >>>
  
  output {
    Array[File] weighted_roc_files = glob("*weighted_roc.tsv.gz")
  }

  runtime {
      docker: "us.gcr.io/broad-dsde-methods/rtg:3.11_091520"
      memory: "10 GB"
      disks: "local-disk " + 30 + " HDD"
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
    touch files.txt
    for roc_file in ~{sep=' ' weighted_roc_files}  ; do
        echo $roc_file >> files.txt
    done
    filename_list=`cat files.txt`
    rtg rocplot ${filename_list} --png=rtgrocplot.png --precision-sensitivity
  >>>

  output {
    File rocplot = "rtgrocplot.png"
  }

  runtime {
      docker: "us.gcr.io/broad-dsde-methods/rtg:3.11_091520"
      memory: "10 GB"
      disks: "local-disk " + 50 + " HDD"
      preemptible: select_first([preemptible_tries, 5])
      cpu: 1
  }
}

