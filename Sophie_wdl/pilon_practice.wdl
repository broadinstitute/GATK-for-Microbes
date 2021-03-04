version 1.0

workflow pilon_practice {
  input {
    File reference
    File bam_file
    File bai_file
    String message
  }

  call pilon_work { input: reference=reference, bam_file=bam_file, bai_file=bai_file, message=message}
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
    pilon --genome ~{reference} --frags ~{bam_file} --variant --outdir . #this puts the files in cromwell-executions/pilon_practice/(most recent run folder)/call-pilon_work/execution/
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
