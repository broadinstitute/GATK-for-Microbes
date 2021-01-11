version 1.0

workflow GenerateRTGData {

  input {
    File truth_fasta
    File baseline_vcf
    Array[File] sample_vcfs
 
    Int? preemptible_tries
  }

  call IndexReference {
    input:
      ref_fasta=truth_fasta
  }

  call BGZipVcf as ProcessBaselineVcf {
    input:
      sample_vcf=baseline_vcf,
      preemptible_tries = preemptible_tries
  }

  scatter (vcf in sample_vcfs) {
      call BGZipVcf {
        input:
          sample_vcf = vcf,
          preemptible_tries = preemptible_tries
    }
    call GATKConcordance {
      input:
        ref_fasta=truth_fasta,
        ref_fasta_index=IndexReference.ref_fasta_fai,
        ref_dict=IndexReference.ref_dict,
        input_vcf=BGZipVcf.vcf_gz,
        input_vcf_idx=BGZipVcf.vcf_gz_index,
        truth_vcf=ProcessBaselineVcf.vcf_gz,
        truth_vcf_idx=ProcessBaselineVcf.vcf_gz_index
    }
 }

  call EvalVcf {
    input:
      truth_fasta=truth_fasta,
      baseline_vcf_gz = ProcessBaselineVcf.vcf_gz,
      baseline_vcf_gz_index = ProcessBaselineVcf.vcf_gz_index,
      sample_vcf_gz_list = BGZipVcf.vcf_gz,
      sample_vcf_gz_index_list = BGZipVcf.vcf_gz_index,
      preemptible_tries = preemptible_tries
}

  call GeneratePlots {
    input:
      weighted_roc_files = EvalVcf.weighted_roc_files,
      preemptible_tries = preemptible_tries
  }

  output {
      Array[File] gatk_concordance_summary_files = GATKConcordance.concordance_tsv
      Array[File] weighted_roc_files = EvalVcf.weighted_roc_files
      Array[File] rtg_summary_txt = EvalVcf.rtg_summary_txt
      File rocplot = GeneratePlots.rocplot
  }
}


task BGZipVcf {
  input {
    File sample_vcf
     # runtime
    Int? preemptible_tries
  }

  Int disk_size = ceil(size(sample_vcf, "GB")) + 20
  String local_vcf = basename(sample_vcf)
  String vcf = basename(local_vcf, ".gz")
  
  command <<<
    ln -s ~{sample_vcf} ~{local_vcf}
    # this will unzip if it is zipped (we specifically need a bgzip, not just gzip)
    gunzip ~{local_vcf} -q
    bgzip ~{vcf}
    tabix -p vcf ~{vcf}.gz
  >>>
  
  output {
    File vcf_gz = "~{vcf}.gz"
    File vcf_gz_index = "~{vcf}.gz.tbi"
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
      ln -s ${sample_basename}_vcfeval/summary.txt ${sample_basename}_summary.txt
    done
    
  >>>
  
  output {
    Array[File] weighted_roc_files = glob("*weighted_roc.tsv.gz")
    Array[File] rtg_summary_txt = glob("*summary.txt")
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

task GATKConcordance {
  input {
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File input_vcf
    File input_vcf_idx
    File truth_vcf
    File truth_vcf_idx
    Int? preemptible_tries
    File? gatk_override
    String? gatk_docker_override
  }

  String output_tsv = basename(basename(input_vcf, ".gz"), ".vcf") + "_concordance_summary.tsv" 
  String base_fasta = basename(ref_fasta)
  String base_fai = basename(ref_fasta_index)
  String base_dict = basename(ref_dict)

  command {
    set -e
    ln -s ~{ref_fasta_index} ~{base_fasta}
    ln -s ~{ref_fasta_index} ~{base_fai}
    ln -s ~{ref_dict} ~{base_dict}

    export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}
    gatk Concordance \
      -R ~{base_fasta} \
      -eval ~{input_vcf} \
      --truth ~{truth_vcf} \
      --summary ~{output_tsv}
  }
  output {
    File concordance_tsv = "~{output_tsv}"
  }
  runtime {
      docker: select_first([gatk_docker_override, "us.gcr.io/broad-gatk/gatk:4.1.7.0"])
      memory: "3 MB"
      disks: "local-disk 20 HDD"
      preemptible: select_first([preemptible_tries, 5])
  } 
}

task IndexReference {
  input {
    File ref_fasta    
    Int? preemptible_tries
  }

  Int disk_size = ceil(size(ref_fasta, "GB") * 2.5) + 20
  String fasta_nopath = basename(ref_fasta)
  String fasta_with_suffix = basename(fasta_nopath, ".gz")
  String basename = basename(basename(fasta_with_suffix, "sta"), ".fa") 
  command <<<
      set -e
      ln -s ~{ref_fasta} ~{fasta_nopath}
      samtools faidx ~{fasta_nopath}
      samtools dict ~{fasta_nopath} > ~{basename}.dict
      ls -al
  >>>
  runtime {
    preemptible: select_first([preemptible_tries, 5])
    memory: "2 GB"
    disks: "local-disk " + disk_size + " HDD"
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.2-1552931386"
  }

  output {
    File ref_dict = "~{basename}.dict"
    File ref_fasta_fai = "~{fasta_with_suffix}.fai"
  }
}

