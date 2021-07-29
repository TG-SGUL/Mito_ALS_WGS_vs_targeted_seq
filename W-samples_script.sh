#!/bin/bash -l
#SBATCH --output=/scratch/users/%u/%j.out
#SBATCH --job-name=Wsamples
#SBATCH --time=0-2:00
#SBATCH --mem=10240
#SBATCH --nodes=1
#SBATCH --ntasks=4

# Abort if there is an error
set -e

# ------------------------
## 1. SETUP
# ------------------------

# 1.a. FILE LOCATIONS
# Use absoloute paths when running as a batch job, and use scratch/ paths when possible 
RAW=/scratch/groups/mito_als_project/WGS_bam_files/batch3/not_processed
FASTQC=/scratch/groups/mito_als_project/WGS_bam_files/batch3/fastqc_and_multiqc

PREALIGN=/scratch/groups/mito_als_project/aligned_data/Exp_22/pre-realignment
ALIGN=/scratch/groups/mito_als_project/aligned_data/Exp_22/alignments

REF_GRCH38=/scratch/groups/mito_als_project/reference/GRCh38-p13.fa

LOGS=/scratch/groups/mito_als_project/documents/experiment_22/logs/W
CONLOG=/scratch/groups/mito_als_project/documents/experiment_22/logs/W/conversion_log.csv

# 1.b. MANUALLY-SET VARIABLES FOR READ GROUP ASSIGNMENT
# What platform were the sequences generated on, in all caps (e.g. ILLUMINA)
PLATFORM="ILLUMINA"
# What model of this platform, in all caps (e.g. MISEQ)
MODEL="HISEQ2000"
# What centre was the sequencing performed at, in all caps (e.g. LIVERPOOL)
CENTRE="ILLUMINA_FASTTRACK"
# What library was used during sequencing, in all caps (e.g )
LIBRARY="TRUSEQ"

# 1.c. MANUALLY-SET VARIABLES FOR TOOLS VERSION LOG
# What version of BWA are you running?
VERSION_BWA="Version 0.7.17-r1188"
# What version of Picard are you running?
VERSION_PICARD="Version:2.25.4"
# What version of VCFlib or VCFfilter are you running?  
VERSION_VCFLIB="vcflib 1.0.2"


# 1.d. CREATE TOOLS VERSION LOG

echo -e "# Versions of tools used\n" > "$LOGS"/W-script_tools_version_log.md
    Vfastqc=$(fastqc --version)
        echo -e "### FastQC\n$Vfastqc\n" >> "$LOGS"/W-script_tools_version_log.md
    Vmultiqc=$(multiqc --version)
        echo -e "### MultiQC\n$Vmultiqc\n" >> "$LOGS"/W-script_tools_version_log.md
    Vbwa="$VERSION_BWA"
        echo -e "### BWA\n$Vbwa\n" >> "$LOGS"/W-script_tools_version_log.md
    Vsamtools=$(samtools version | head -n 1)
        echo -e "### Samtools\n$Vsamtools\n" >> "$LOGS"/W-script_tools_version_log.md
    Vpicard="$VERSION_PICARD"
        echo -e "### Picard CollectInsertSizeMetrics\n$Vpicard\n" >> "$LOGS"/W-script_tools_version_log.md  
    Vbedtools=$(bedtools --version)
        echo -e "### BEDtools\n$Vbedtools\n"  >> "$LOGS"/W-script_tools_version_log.md
    Vfreebayes=$(freebayes --version)
        echo -e "### Freebayes\n$Vfreebayes\n"  >> "$LOGS"/W-script_tools_version_log.md  
    Vbcftools=$(bcftools --version | head -n 2)
        echo -e "### BCFtools\n$Vbcftools\n"  >> "$LOGS"/W-script_tools_version_log.md 
    Vvcffilter="$VERSION_VCFLIB"
        echo -e "### VCFfilter\n$Vvcffilter\n"  >> "$LOGS"/W-script_tools_version_log.md   
    Vgatk=$(gatk --version)
        echo -e "### GATK\n$Vgatk\n"  >> "$LOGS"/W-script_tools_version_log.md

# 1.e. CREATE READ-COUNT LOG
echo "Sample_ID, chrM_reads_aligning_whole_genome, rCRS_aligning_reads" > "$LOGS"/W_read_counts.csv
read_count_log="$LOGS"/W_read_counts.csv

# 1.f. CREATE BAM TO FASTA CONVERSION LOG
echo "ID,paired_read_count,unpaired_read_count,fastq_R1_reads,fastq_R2_reads" > "$CONLOG"

# ------------------------
## 2. CONVERSION TO FASTQ
# ------------------------

for f in "$RAW"/*.bam
        do

        # Generate sample name (ID)
        ID=W_$(basename "$f" | sed -e 's/\.bwa\.sorted\.marked\.bam$//')
        
        # Generate a sample folder for each sample
        # Note that this will cause an ERROR if the sample folder already exists!
        mkdir "$PREALIGN"/"$ID"
        SAMPLEFOLDER="$PREALIGN"/"$ID"

        echo "$ID applying filtering"

        # Filter to only chrM reads 
        samtools view -h \
        -O sam \
        -o "$SAMPLEFOLDER"/"$ID"_chrM.sam \
        "$f" \
        chrM

        # Inclusionary flag for paired and mapped reads(-f 2)
        # Output with header (-h), as sam file (-O)
        samtools view -h -f 2 \
        -O sam \
        -o "$SAMPLEFOLDER"/"$ID"_pre-realignment_paired.sam \
        -U "$SAMPLEFOLDER"/"$ID"_pre-realignment_unpaired.sam \
        "$SAMPLEFOLDER"/"$ID"_chrM.sam

        # Convert to fastq files for R1 and R2
        echo "$ID converting to fastq"
        
        # BEDTOOLS METHOD
        samtools sort -n \
        -O bam \
        -o "$SAMPLEFOLDER"/"$ID"_pre-realignment_sorted_paired.bam \
        "$SAMPLEFOLDER"/"$ID"_pre-realignment_paired.sam
        
        bedtools bamtofastq -i "$SAMPLEFOLDER"/"$ID"_pre-realignment_sorted_paired.bam \
        -fq "$SAMPLEFOLDER"/"$ID"_R1.fastq \
        -fq2 "$SAMPLEFOLDER"/"$ID"_R2.fastq

        echo "$ID updating conversion log"
        # Write information to the conversion log
        # Number of reads in R1 and R2
        R1LINES=$(wc -l "$SAMPLEFOLDER"/"$ID"_R1.fastq |cut -d " " -f 1)
        R1COUNT=$(($R1LINES / 4))
        R2LINES=$(wc -l "$SAMPLEFOLDER"/"$ID"_R2.fastq |cut -d " " -f 1)
        R2COUNT=$(($R2LINES / 4))
        # Paired and unpaired reads
        PAIRED=$(samtools view -c "$SAMPLEFOLDER"/"$ID"_pre-realignment_paired.sam)
        UNPAIRED=$(samtools view -c "$SAMPLEFOLDER"/"$ID"_pre-realignment_unpaired.sam)
        echo "$ID,$PAIRED,$UNPAIRED,$R1COUNT,$R2COUNT" >> "$CONLOG"

        echo "$ID cleanup"
        # Clean up: remove pre-realignment sam and zip files
        rm "$SAMPLEFOLDER"/"$ID"_pre-realignment_paired.sam
        rm "$SAMPLEFOLDER"/"$ID"_chrM.sam

    done

#------------------------------------
# 3. FASTQ ASSESSMENT OF RAW FILES
#------------------------------------

# 3.a. GENERATE FASTQ FILES 
for f in "$PREALIGN"/*/*_R1.fastq
    do
    # Generate R1 unzipped filename
    R1=$(echo $f)
    # Generate R2 filenames - zipped and unzipped
    R2=$(echo $f | sed 's/_R1/_R2/' )

    # Run fastqc
    fastqc -t 6 -o "$FASTQC" "$f"
    fastqc -t 6 -o "$FASTQC" "$R2"
done

# 3.b. GENERATE MULTIQC FILE
# Enter folder with fastqc files, execute MultiQC, return to start location
pushd "$FASTQC"
multiqc .
popd

echo "
    ------------------------------------------------------------------
    Quality reports generated. 
    Examine MultiQC report and decide whether to move on to alignment.
    ------------------------------------------------------------------"
#

# ------------------------
## 4. ALIGNMENT
# ------------------------

for f in "$PREALIGN"/*/*_R1.fastq
    do
    # Extract ID
    ID=$(basename "$f" | sed -e 's/_R1.fastq$//')
    # Generate R2 filename
    R2=$(echo $f | sed 's/_R1/_R2/')
    
    # Make each sample their own data folder
    # THIS WILL INDUCE AN ERROR IF THIS FOLDER ALREADY EXISTS
    mkdir "$ALIGN"/"$ID"
    SAMPLEFOLDER="$ALIGN"/"$ID"

    # Extract read group information
    echo "$ID - gathering @RG information"
    ID2="$ID".$(head -n 1 "$f" | cut -d : --output-delimiter=. -f 3,4)
    SM="$ID"
    LN=$(head -n 1 "$f" | cut -d : --output-delimiter=. -f 3,4)
    PU=$(head -n 1 "$f" | cut -d : --output-delimiter=. -f 3,4)
    # Additional read group information for targeted sequences
    PL=$PLATFORM
    PM=$MODEL
    CN=$CENTRE
    LB=$LIBRARY

    # Perform Burrows-Wheeler alignment to GRCh38
    echo "$ID performing alignment"
    bwa mem \
    -v 2 \
    -t "$(nproc)" \
    -R "@RG\tID:${ID2}\tSM:${SM}\tPL:${PL}\tPM:${PM}\tCN:${CN}\tLB:${LB}\tLN:${LN}\tPU:${PU}" \
    $REF_GRCH38 \
    "$f" "$R2" > "$SAMPLEFOLDER"/"$ID"_aligned-GRCh38.sam
    
    # Convert to bam format
    samtools view -h -b -o "$SAMPLEFOLDER"/"$ID"_aligned-GRCh38_with_discordant.bam \
    "$SAMPLEFOLDER"/"$ID"_aligned-GRCh38.sam
    
    # Remove discordant reads immediately after alignment, as additional filtering can disrupt this flag
    samtools view -b -f 2 -o "$SAMPLEFOLDER"/"$ID"_aligned-GRCh38.bam \
    "$SAMPLEFOLDER"/"$ID"_aligned-GRCh38_with_discordant.bam
    
    # Sort and index aligned file
    samtools sort -o "$SAMPLEFOLDER"/"$ID"_aligned-GRCh38_sorted.bam \
    -O bam \
    -@ "$(nproc)" \
    "$SAMPLEFOLDER"/"$ID"_aligned-GRCh38.bam 
    # Index prefilter files
    samtools index -@ "$(nproc)" \
    "$SAMPLEFOLDER"/"$ID"_aligned-GRCh38_sorted.bam 

    # Extract reads which match chrM. 
    # NC_012920.1 is the NCBI nomenclature for chrM.
    echo "$ID - alignment finished, extracting reads from chrM"
    whole_genome_readcount=$(samtools view -c "$SAMPLEFOLDER"/"$ID"_aligned-GRCh38_sorted.bam) 

    # Extracting reads
    samtools view -h -o "$SAMPLEFOLDER"/"$ID"_rCRS.sam \
    "$SAMPLEFOLDER"/"$ID"_aligned-GRCh38_sorted.bam NC_012920.1
    
    # Converting to bam format
    samtools view -h -b -o "$SAMPLEFOLDER"/"$ID"_rCRS.bam \
    "$SAMPLEFOLDER"/"$ID"_rCRS.sam

    # Output info to read count log
    rCRS_only_readcount=$(samtools view -c "$SAMPLEFOLDER"/"$ID"_rCRS.bam)
    echo "$ID, $whole_genome_readcount, $rCRS_only_readcount" >> "$read_count_log"
    
    # Convert chrM to bam, sort and index
    samtools sort -o "$SAMPLEFOLDER"/"$ID"_rCRS_sorted.bam -O bam -@ "$(nproc)" "$SAMPLEFOLDER"/"$ID"_rCRS.bam
    samtools index -@ "$(nproc)" "$SAMPLEFOLDER"/"$ID"_rCRS_sorted.bam

    # Clean-up the larger .sam files and intermediate bam and bai files
    rm "$SAMPLEFOLDER"/"$ID"_rCRS.sam
    rm "$SAMPLEFOLDER"/"$ID"_rCRS.bam
    rm "$SAMPLEFOLDER"/"$ID"_aligned-GRCh38.bam
    rm "$SAMPLEFOLDER"/"$ID"_aligned-GRCh38.sam
    rm "$SAMPLEFOLDER"/"$ID"_aligned-GRCh38_sorted.bam.bai
    rm "$SAMPLEFOLDER"/"$ID"_aligned-GRCh38_with_discordant.bam

done

echo "
    -------------------------------------------------------------------------------
    W-sample alignment complete
    To move on to filtering, QC, and variant calling, run the All-samples_script.sh
    --------------------------------------------------------------------------------"