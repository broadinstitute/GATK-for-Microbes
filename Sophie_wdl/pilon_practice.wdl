version 1.0

workflow pilon_practice {
  input {
    File reference
    File bam_file
    #File bai_file
    String message
  }

  String base = basename(bam_file)
  File bai_file = base + ".bai"

  call pilon_work { input: reference=reference, bai_file=bai_file, bam_file=bam_file, message=message}
  call filtering { input: vcf=pilon_work.vcf}
}

task pilon_work {
  input {
    File reference
    File bam_file
    File bai_file
    String message
  }
  command {
    echo '${message}!'
    pilon -Xmx8G --genome ~{reference} --frags ~{bam_file} --variant --tracks --outdir . #this puts the files in cromwell-executions/pilon_practice/(most recent run folder)/call-pilon_work/execution/
  }
  output {
    #Pilon otput files
    File vcf = "pilon.vcf"
    File fasta = "pilon.fasta"
  }
  runtime {
   docker: 'quay.io/biocontainers/pilon:1.23--2'
  }
}

task filtering {
  input {
    File vcf
  }
  command {
    bcftools filter -e'ALT="."' --output filtered.vcf --output-type v ~{vcf}
  }
  output {
    #Pilon otput files
    File vcf_filtered = "filtered.vcf"
  }
}


