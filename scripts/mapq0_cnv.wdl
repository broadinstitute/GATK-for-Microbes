version 1.0

workflow mapq0_cnv {
  input {
    File bam_file
    Int? distance = 0
    Int? cutoff = 50
  }

  String base = basename(bam_file)
  File bai_file = base + ".bai"

  call mapq0_cnv_calls { input: bai_file=bai_file, bam_file=bam_file, distance=distance, cutoff=cutoff}
}

task mapq0_cnv_calls {
  input {
    File bam_file
    File bai_file
    Int distance
    Int cutoff
  }
  command {
    python mapq0_cnv.py -i ~{bam_file} -d ~{distance} -c ~{cutoff} -o 'output.bed'
  }
  output {
    File output = 'output.bed'
  }
  runtime {
   #docker: 'quay.io/biocontainers/pilon:1.23--2'
  }
}
