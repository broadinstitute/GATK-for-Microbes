version 1.0

workflow convertSamToFastq {

  input {
	File inputBam
	String sampleName
	Int memoryGb
	Int diskSpaceGb
  }
    
	call samToFastq {
	  input:
		inputBam = inputBam,
		sampleName = sampleName,
		memoryGb = memoryGb,
		diskSpaceGb = diskSpaceGb
	}
}

task samToFastq {
	input {
		File inputBam
		String sampleName
		Int memoryGb
		Int diskSpaceGb
	}

	command <<<
		java -jar /usr/gitc/picard.jar SamToFastq \
		TMP_DIR=. \
		INPUT=~{inputBam} \
		FASTQ=~{sampleName}.1.fastq.gz \
		INTERLEAVE=false \
		SECOND_END_FASTQ=~{sampleName}.2.fastq.gz \
		INCLUDE_NON_PF_READS=true \
		CLIPPING_ATTRIBUTE=XT \
		CLIPPING_ACTION=2 \
		UNPAIRED_FASTQ=~{sampleName}.unpaired.fastq.gz
	>>>

	output {
		File firstEndFastq = "~{sampleName}.1.fastq.gz"
		File secondEndFastq = "~{sampleName}.2.fastq.gz"
		File unpairedFastq = "~{sampleName}.unpaired.fastq.gz"
	}

	runtime {
		docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1500064817"
		memory: "~{memoryGb} GB"
		cpu: "1"
		disks: "local-disk ~{diskSpaceGb} HDD"
		
	}
}
