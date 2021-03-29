version 1.0

workflow mummer_practice {
  input {
    File reference
    File alignment
  }

  call all { input: reference=reference, alignment=alignment}
}

task all {
  input {
    File reference
    File alignment
  }
  command {
    nucmer --mum -p alignment ~{reference} ~{alignment}
    delta-filter -1 -i 95 -l 100 alignment.delta > alignment.fil.delta
    show-coords alignment.fil.delta > alignment.fil.coords
    show-snps -C alignment.fil.delta > alignment.fil.delta.snps
  }
  output {
    File delta = "alignment.delta"
    File delta_filtered = "alignment.fil.delta"
    File filtered_cords = "alignment.fil.coords"
    File filtered_snps = "alignment.fil.delta.snps"
  }
  runtime {
   docker: 'quay.io/biocontainers/mummer:3.23--pl526he1b5a44_11'
  }
}