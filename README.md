# GATK-for-Microbes


cromshell submit MicrobialGenomePipeline.wdl align-inputs.json

### GATK Best Practices for Variant Calling in Microbial Genomes

A reproducible pipeline for SNP and Indel variant calling in microbial whole-genome sequencing data.  

This pipeline is in the **alpha stage** of development. We have done our intial benchmarking and have seen promising results with sensitivity and precision. We are  sharing this pipeline with the community to get feedback and make further improvements where required. 

### Summary 

Genome Analysis Toolkit (GATK) is the industry standard for human variant discovery and genotyping. Due to its highly effective variant discovery methods, researchers are adopting GATK for microbes even though microbial genomes are radically different from human genomes and GATK has not been optimized for it. To optimize GATK for microbes and provide researchers with robust resources to study and combat life threatening pathogens, we developed a scalable, portable, and reproducible workflow that calls high-quality filtered variants using short reads sequencing data. We optimized for low allele frequencies, varying read depths, and sequencing and mapping errors typical of bacterial data, resulting in improved sensitivity and precision. For improved coverage across circular bacterial genomes, we developed a tool that calls variants on reads spanning the breakpoint at which they are linearized. We started these efforts with bacteria and will expand it to cover other microbial pathogens.

## Workflows
![pipeline](https://drive.google.com/uc?export=view&id=12RMqb-jkw6RDGgEhS1JoEr7kwQ5ycd5M)


### 1. MicrobialGenomePipeline

**What does it do?**     
This WDL workflow performs SNP/Indel variant calling on WGS input data.

**What data does it require as input?**  
This workflow accepts WGS FASTQ or BAM files as inputs and closest reference FASTA files. 

**What does it output?** .    
A filtered VCF file of SNP/Indel calls.         

**Sample data description and location**    
An example bam input of *Mycobacterium tuberculosis* is provided in the workspace data model for testing.    

**Reference data description and location**  
The required and optional references and resources for the Tools are included in the Workspace Data table.       

**Time and cost estimates**    

| Participant | Size | Time | Cost $ |
| :------------------: | :----------------: | :------: | :--------: |
| mtb | 675.12 MB | 0:16:00 | <0.01 |

**Software Version Notes**   
GATK 4.1.9.0  

### 2. MicrobialConcordancePipeline

**What does it do?**     
This WDL workflow performs concordance on one or more variant files (vcf) against a truth vcf.

**What data does it require as input?**  
This workflow accepts one or more samples and truth vcf files as inputs and closest reference FASTA files. 

**What does it output?** .    
Results of concordance comparison, summary metrics, and ROC curve data files.      

Cost and Time will vary, view [Controlling-Cloud-costs-sample-use-cases](https://support.terra.bio/hc/en-us/articles/360029772212) for further details.  
Users can also use [Google's BigQuery](https://software.broadinstitute.org/firecloud/documentation/article?id=11788) for task level calculation. 

### 3. DisplayConcordance Jupyter Notebook

We have also provided a notebook which takes as input the results of the MicrobialConcordancePipeline and plots the concordance comparison as shown in the image below. This notebook helps us streamline and quickly test new improvements we make to the pipeline against a) previous versions of the pipeline and/or b) other tools/pipelines.

![notebook](https://drive.google.com/uc?export=view&id=1Z_39Gv6LbvDqa1obqoTfR7PIYoKiqhFC)

---

### Contact Information  
For questions about this workspace please reach out to Bhanu Gandham at bgandham@broadinstitute.org

This material is provided by the GATK Team. Please post any questions or concerns regarding the GATK tool to the [GATK forum](https://gatk.broadinstitute.org/hc/en-us/community/topics)

### License  
**Copyright Broad Institute, 2019 | BSD-3**  
All code provided in this workspace is released under the WDL open source code license (BSD-3) (full license text at https://github.com/openwdl/wdl/blob/master/LICENSE). Note however that the programs called by the scripts may be subject to different licenses. Users are responsible for checking that they are authorized to run all programs before running these tools.

### Workspace Change Log
| Date | Change | Author |
| --- | --- | --- |
|  2020-11-12 | Initial Setup | Bhanu Gandham and Andrea Haessly |




