import sys, os


### This should be implemented with argparse to make it prettier
### But for nowI'll just set these via sys.argv

##################
### Initialize ###
##################

if len(sys.argv) < 5:
	print "Usage: python Pipeline_Master.py <Sample_ID> <Working_Directory> <R1.fastq> <R2.fastq>"
	print "Use full filepaths for the working directory and fastq files"
	sys.exit()
sampleID = sys.argv[1]
workingDir = sys.argv[2]
R1fastq = sys.argv[3]
R2fastq = sys.argv[4]

print "The working directory is: %s"%workingDir
print "The sample ID is: %s"%sampleID
print "The fastq files you're working with are: %s\n%s"%(R1fastq,R2fastq)

#maybe include an option for processors and genome version?
numProcessors = 32
Memory = "100G"
runTime = "4:0:0"

#######################
### Write the Script###
#######################

shellScriptFile = open('%s_FullPipeline.sh'%sampleID,'w')
shellScriptFile.write('#!/bin/bash\n')



#Job name
shellScriptFile.write('#$ -N PrimaryPipeline_%s\n'%sampleID)
#Export environment variables
shellScriptFile.write('#$ -V\n')
#set location for log files
shellScriptFile.write('#$ -o %s%s.e\n'%(workingDir,sampleID))
#email on job abort
shellScriptFile.write('#$ -m a\n#$ -M prichmond@cmmt.ubc.ca\n')
#SGE resource requirements
shellScriptFile.write('#$ -l h_rt=%s,h_vmem=%s\n'%(runTime,Memory))
#Parallel processing
shellScriptFile.write('#$ -pe smp %d\n'%numProcessors)
shellScriptFile.write('\n\nexport PARALLEL=$NSLOTS\nexport OMP_NUM_THREADS=$NSLOTS\n')

#Set the variables for working directory, sampleID, and fastq full filepath
shellScriptFile.write("\nSAMPLE_ID=\'%s\'\n"%sampleID)
shellScriptFile.write("WORKING_DIR=\'%s\'\n"%workingDir)
shellScriptFile.write("FASTQR1=\'%s\'\n"%R1fastq)
shellScriptFile.write("FASTQR2=\'%s\'\n"%R2fastq)
shellScriptFile.write("BOWTIE2_INDEX=\'/mnt/data/GENOMES/hg19/hg19\'\n")


#Map with Bowtie2
shellScriptFile.write("\n#Map with Bowtie2\n")
shellScriptFile.write("/opt/tools/bowtie2-2.2.6/bowtie2  -x $BOWTIE2_INDEX -1 $FASTQR1 -2 $FASTQR2 -S $WORKING_DIR$SAMPLE_ID\'_bowtie2.sam\'  -p $NSLOTS --very-sensitive --score-min L,-10,-0.3 -X 1000 --met-stderr --rg-id $SAMPLE_ID --rg \"SM:$SAMPLE_ID\" 2> $WORKING_DIR$SAMPLE_ID\'_bowtie2.stderr\'\n\n")

#Convert to binary, sort, and index
shellScriptFile.write("#Convert to binary, sort, and index\n")
shellScriptFile.write("/opt/tools/samtools/samtools view -bS $WORKING_DIR$SAMPLE_ID\'_bowtie2.sam\' | /opt/tools/samtools/samtools sort -m100000000000 - $WORKING_DIR$SAMPLE_ID\'_bowtie2.sorted\'\n")
shellScriptFile.write("/opt/tools/samtools/samtools index $WORKING_DIR$SAMPLE_ID\'_bowtie2.sorted.bam\'\n")

#Remove duplicates
shellScriptFile.write("\n#Remove Duplicates\n")
shellScriptFile.write("/opt/tools/jdk1.7.0_79/bin/java -jar /opt/tools/picard-tools-1.139/picard.jar MarkDuplicates I=$WORKING_DIR$SAMPLE_ID\'_bowtie2.sorted.bam\' O=$WORKING_DIR$SAMPLE_ID\'_bowtie2_dupremoved.sorted.bam\' REMOVE_DUPLICATES=true M=$WORKING_DIR$SAMPLE_ID\'_bowtie2_DuplicateResults.txt\'\n")
shellScriptFile.write("/opt/tools/samtools/samtools index $WORKING_DIR$SAMPLE_ID\'_bowtie2_dupremoved.sorted.bam\'\n")

#Realign
shellScriptFile.write("\n#Realign\n")
shellScriptFile.write("/opt/tools/jdk1.7.0_79/bin/java -jar /opt/tools/GATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R /mnt/data/GENOMES/hg19/FASTA/hg19.fa -minReads 5 -I $WORKING_DIR$SAMPLE_ID\'_bowtie2_dupremoved.sorted.bam\' -o $WORKING_DIR$SAMPLE_ID\'_bowtie2_indelsites.intervals\' -nt $NSLOTS\n")
shellScriptFile.write("/opt/tools/jdk1.7.0_79/bin/java -jar /opt/tools/GATK/GenomeAnalysisTK.jar -T IndelRealigner -model USE_READS -R /mnt/data/GENOMES/hg19/FASTA/hg19.fa -targetIntervals $WORKING_DIR$SAMPLE_ID\'_bowtie2_indelsites.intervals\' -I $WORKING_DIR$SAMPLE_ID\'_bowtie2_dupremoved.sorted.bam\' -o $WORKING_DIR$SAMPLE_ID\'_bowtie2_dupremoved_realigned.sorted.bam\'\n")
shellScriptFile.write("/opt/tools/samtools/samtools index $WORKING_DIR$SAMPLE_ID\'_bowtie2_dupremoved_realigned.sorted.bam\'\n")

#SNP calling
shellScriptFile.write("\n#SNP Calling\n")
shellScriptFile.write("/opt/tools/jdk1.7.0_79/bin/java -jar /opt/tools/GATK/GenomeAnalysisTK.jar -T UnifiedGenotyper -nt $NSLOTS -glm BOTH -R /mnt/data/GENOMES/hg19/FASTA/hg19.fa -I $WORKING_DIR$SAMPLE_ID\'_bowtie2_dupremoved_realigned.sorted.bam\' -o $WORKING_DIR$SAMPLE_ID\'_bowtie2_dupremoved_realigned_UnifiedGenotyper.vcf\'\n")

