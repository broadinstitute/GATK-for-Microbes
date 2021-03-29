version 1.0

workflow minimap {
  input {
    File reference
    File alignment
  }

  call minimap_work { input: reference=reference, alignment=alignment}
}

task minimap_work {
  input {
    File reference
    File alignment
  }
  command {
    minimap2 -cx asm5 --cs ~{reference} ~{alignment} | sort -k6,6 -k8,8n | paftools.js call -f ~{reference} - > out.vcf
  }
  output {
    File vcf = "out.vcf"
  }
}
