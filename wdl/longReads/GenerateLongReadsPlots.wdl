version 1.0

workflow GenerateLongReadsPlots {
	input {
        File reads
    }

        call GenerateRPlots {
        input: 
            reads = reads
    }

    output {
        File zmw_frequency_plot = GenerateRPlots.zmw_frequency_plot
    }
}    

    task GenerateRPlots{
    input {
        File reads
    }
    command <<<
  set -e
  zcat ~{reads} | zgrep '@' | awk -F '/' '{ print $2 }' > zmw_list

  R --vanilla <<RSCRIPT
  library('tidyverse')
  library('dplyr')
  zmw = read.table('zmw_list')
  zmw_count = zmw %>% group_by(V1) %>% count(V1)
  count = unlist(zmw_count[,'n'])
  head(count)
  pdf("zmw_frequency_plot.pdf")
  hist(count, xlim = c(0,50), breaks=500)
  dev.off()

  RSCRIPT
  >>>

  output {
        File zmw_frequency_plot = "zmw_frequency_plot.pdf"
    }
    runtime {
      docker: "rocker/tidyverse:4.0.3"
      memory: "2 GB"
      disks: "local-disk " + 10 + " HDD"
      preemptible: 3
      cpu: 1
  }
}
