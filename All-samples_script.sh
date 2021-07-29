#!/bin/bash -l
#SBATCH --output=/scratch/users/%u/%j.out
#SBATCH --job-name=AllSamples
#SBATCH --time=0-2:00
#SBATCH --mem=16384
#SBATCH --nodes=1
#SBATCH --ntasks=4

# Abort if there is an error
set -e

# ------------------------
## 1. SETUP
# ------------------------

# 1.a. FILE LOCATIONS
# Use absoloute paths when running as a batch job, and use scratch/ paths when possible 

ALIGN=/scratch/groups/mito_als_project/aligned_data/Exp_22/alignments
RESULTS=/scratch/groups/mito_als_project/aligned_data/Exp_22/metrics

VARS=/scratch/groups/mito_als_project/variant_calls/exp_22/vcf_files
VARTABLE=/scratch/groups/mito_als_project/variant_calls/exp_22/variant_tables

REF_RCRS=/scratch/groups/mito_als_project/reference/rCRS.fa

# # ------------------------
# ## 5. ALIGNMENT FILTERING 
# # ------------------------

# # 5.a. GENERATE PRE-FILTER METRICS

# # Create summary files
#     prefilter_flagstat="$RESULTS"/prefilter_flagstat.csv
#     echo "ID,reads,secondary,supplementary,duplicates,mapped,paired_in_sequencing,read1,read2,properly_paired,itself_and_mate_mapped,singletons,mate_mapped_elsewhere,mate_mapped_elsewhere_mapq>=5" > $prefilter_flagstat

#     prefilter_coverage="$RESULTS"/prefilter_coverage.csv
#     echo "ID,rname,startpos,endpos,numreads,covabses,coverage,meandepth,meanbaseq,meanmapq" > $prefilter_coverage
    
#     prefilter_insertsize="$RESULTS"/prefilter_insertsize.csv
#     echo "ID,median_insertsize,mode_insertsize,median_absoloute_deviation,min_insertsize,max_insertsize,mean_insertsize,stdev,read_pairs,pair_orientation,width_of_10_percent,width_of_20_percent,width_of_30_percent,width_of_40_percent,width_of_50_percent,width_of_60_percent,width_of_70_percent,width_of_80_percent,width_of_90_percent,width_of_95_percent,width_of_99_percent,SAMPLE,LIBRARY,READ_GROUP" > $prefilter_insertsize

#     prefilter_summary="$RESULTS"/prefilter_summary.csv

# # Generate pre-filtering summary stats
#     for f in "$ALIGN"/*/*_rCRS_sorted.bam
#         do
#         # create filename ID without extension
#         ID=$(basename -s _rCRS_sorted.bam $f)
#         # Assign sample folder name
#         SAMPLEFOLDER="$ALIGN"/"$ID"

#         # perform flagstats and output to prefilter_flagstat
#         echo -e "\n$ID flagstats ongoing"
#         stats=$(samtools flagstat "$f" | tr '\n' ','| sed "s/,$/\n/")
#         echo "$ID,$stats" >> $prefilter_flagstat

#         # generate coverage statistics
#         echo "$ID calculating coverage"
#         coverage=$(samtools coverage -H -r NC_012920.1 "$f" | tr '	' ',')
#         echo "$ID,$coverage" >> $prefilter_coverage

#         # generate insert size metrics
#         echo "$ID collecting insert sizes"
#         # The varaible inserts is one line. CollectInsertSizeMetrics generates a long text file - pipe to grep to extract the 2 relevant lines - pipe to sed to replace spaces with commas. Append this csv-formatted line to the summary file. 
#         picard CollectInsertSizeMetrics -I "$f" \
#         -O "$RESULTS"/insertsize/"$ID"_insertsize.txt \
#         -H "$RESULTS"/insertsize/"$ID"_insertsize.pdf \
#         --INCLUDE_DUPLICATES true \
#         -QUIET true -VERBOSITY ERROR
#         inserts=$(grep -A 1 MEDIAN_INSERT_SIZE "$RESULTS"/insertsize/"$ID"_insertsize.txt| tail -n +2 | sed 's/\s/,/g')
#         echo "$ID,$inserts" >> "$prefilter_insertsize"

#     done

#     # Generate pre-filtering summary file
#     echo -e "\n--------------------------------------\nGenerating pre-filtering summary file\n--------------------------------------\n"
    
#     # Create concise flagstats summary
#         # Remove the header (using tail) as it will no longer align to the data if kept. Pipe to sed - replace spaces with commas. Pipe to cut - take relevant lines into a temporary file.
#         tail -n +2 "$RESULTS"/prefilter_flagstat.csv | sed 's/\s/,/g' | cut -d ',' -f 1,2,12,16,20,24,31,37,41,45,53,61,68,78 > "$RESULTS"/ftemp.csv
#         # Recreate headers
#         HEADER="ID,reads,secondary,supplementary,duplicates,mapped,paired_in_sequencing,read1,read2,properly_paired,itself_and_mate_mapped,singletons,mate_mapped_elsewhere,mate_mapped_elsewhere_mapq>=5"
#         # Compile into concise file
#         echo "$HEADER" > "$RESULTS"/prefilter_flagstat_concise.csv
#         cat "$RESULTS"/ftemp.csv >> "$RESULTS"/prefilter_flagstat_concise.csv
#         # Remove the temporary file and the verbose flagstat file
#         rm "$RESULTS"/ftemp.csv
#         rm "$RESULTS"/prefilter_flagstat.csv

#     # Join pre-QC stats files into one summary file
#     join -t , --header "$RESULTS"/prefilter_flagstat_concise.csv $prefilter_coverage > "$RESULTS"/temp_summary.csv
#     join -t , --header "$RESULTS"/temp_summary.csv $prefilter_insertsize > $prefilter_summary
#     rm "$RESULTS"/temp_summary.csv
# #

# # 5.b. FILTER ALIGNMENTS
#     for f in "$ALIGN"/*/*_rCRS_sorted.bam
#         do

#         # create filename ID without extension
#         ID=$(basename -s _rCRS_sorted.bam $f)
#         # Assign sample folder name
#         SAMPLEFOLDER="$ALIGN"/"$ID"

#         echo "$ID filtering"
         
#         # Filter, using a mapq baseline of 20 and an exclusionary bit flag of 2052 (unmapped and/or supplementary reads)
#         samtools view -b -q 20 -F 2052 \
#             -o "$SAMPLEFOLDER"/"$ID"_rCRS_sorted_filtered.bam \
#             -U "$SAMPLEFOLDER"/"$ID"_rCRS_filtered_out.bam \
#             -@ "$(nproc)" \
#             "$f"
#         # Index the sorted_filtered file
#         samtools index "$SAMPLEFOLDER"/"$ID"_rCRS_sorted_filtered.bam

#         # Clean-up: remove intermediate files
#         rm "$SAMPLEFOLDER"/"$ID"_rCRS_sorted.bam.bai
#         # Clean-up: compress files that are retained
#         bgzip -f "$SAMPLEFOLDER"/"$ID"_rCRS_sorted.bam
#         bgzip -f "$SAMPLEFOLDER"/"$ID"_rCRS_filtered_out.bam
#     done

# #

# # 5.c. GENERATE POST-FILTER METRICS
# # Generate summary files
#     postfilter_flagstat="$RESULTS"/postfilter_flagstat.csv
#     echo "ID,reads,secondary,supplementary,duplicates,mapped,paired_in_sequencing,read1,read2,properly_paired,itself_and_mate_mapped,singletons,mate_mapped_elsewhere,mate_mapped_elsewhere_mapq>5" > $postfilter_flagstat

#     postfilter_coverage="$RESULTS"/postfilter_coverage.csv
#     echo "ID,rname,startpos,endpos,numreads,covabses,coverage,meandepth,meanbaseq,meanmapq" > $postfilter_coverage

#     postfilter_insertsize="$RESULTS"/postfilter_insertsize.csv
#     echo "ID,median_insertsize,mode_insertsize,median_absoloute_deviation,min_insertsize,max_insertsize,mean_insertsize,stdev,read_pairs,pair_orientation,width_of_10_percent,width_of_20_percent,width_of_30_percent,width_of_40_percent,width_of_50_percent,width_of_60_percent,width_of_70_percent,width_of_80_percent,width_of_90_percent,width_of_95_percent,width_of_99_percent,SAMPLE,LIBRARY,READ_GROUP" > $postfilter_insertsize

#     postfilter_summary="$RESULTS"/postfilter_summary.csv

# # Generate post-filter summary stats
#     for f in "$ALIGN"/*/*_rCRS_sorted_filtered.bam
#         do
#         # create filename ID without extension
#         ID=$(basename -s _rCRS_sorted_filtered.bam $f)
#         # Assign sample folder name
#         SAMPLEFOLDER="$ALIGN"/"$ID"

#         # perform flagstats, remformat the data to one line, output postfilter_flagstat
#         echo -e "/n$ID flagstats ongoing" 
#         stats=$(samtools flagstat "$f" | tr '\n' ','| sed "s/,$/\n/")
#         echo "$ID,$stats" >> $postfilter_flagstat

#         # generate coverage statistics for the rCRS, replacing any spaces with commas
#         echo "$ID calculating coverage" 
#         coverage=$(samtools coverage -r NC_012920.1 -H "$f" | tr '	' ',')
#         echo "$ID,$coverage" >> $postfilter_coverage
        
#         # generate insert size metrics
#         echo "$ID collecting insert sizes"
#         picard CollectInsertSizeMetrics -I "$f" \
#         -O "$RESULTS"/insertsize/"$ID"_post_insertsize.txt \
#         -H "$RESULTS"/insertsize/"$ID"_post_insertsize.pdf \
#         --INCLUDE_DUPLICATES true \
#         -QUIET true -VERBOSITY ERROR
#         # generate one-line csv-file compatible summary
#         inserts=$(grep -A 1 MEDIAN_INSERT_SIZE "$RESULTS"/insertsize/"$ID"_post_insertsize.txt| tail -n +2 | sed 's/\s/,/g')
#         echo "$ID,$inserts" >> "$postfilter_insertsize"
    
#     done

# # Generate summary file
#     echo -e "\n--------------------------------------\nGenerating post-filtering summary file\n--------------------------------------\n"
    
#     # Create concise flagstats summary
#     # Remove the header as it will no longer align to the data if kept (using tail). Pipe to sed - replace spaces with commas. Pipe to cut - take relevant lines into a temporary file.
#     # This file has one fewer colunm than the pre-qc equivalent. For an explanation see the associated .md document.
#     tail -n +2 "$RESULTS"/postfilter_flagstat.csv | sed 's/\s/,/g' | cut -d ',' -f 1,2,12,16,20,24,31,37,41,45,53,61,68 > "$RESULTS"/temp.csv
#     # Recreate headers, as these would no longer align to the data if kept
#     HEADER="ID,reads,secondary,supplementary,duplicates,mapped,paired_in_sequencing,read1,read2,properly_paired,itself_and_mate_mapped,singletons,mate_mapped_elsewhere_mapq>20"
#     # Compile into concise file
#     echo "$HEADER" > "$RESULTS"/postfilter_flagstat_concise.csv
#     cat "$RESULTS"/temp.csv >> "$RESULTS"/postfilter_flagstat_concise.csv
#     # Remove the temporary file
#     rm "$RESULTS"/temp.csv
#     # Remove verbose flagstat file
#     rm "$RESULTS"/postfilter_flagstat.csv

#     # Join post-QC stats files into one summary file
#     join -t , --header "$RESULTS"/postfilter_flagstat_concise.csv $postfilter_coverage > "$RESULTS"/temp_summary.csv
#     join -t , --header "$RESULTS"/temp_summary.csv $postfilter_insertsize > $postfilter_summary
#     rm "$RESULTS"/temp_summary.csv
# #

# # Generate coverage files
# for f in "$ALIGN"/*/*_rCRS_sorted_filtered.bam
#     do
    
#     # create filename ID without extension
#     ID=$(basename -s _rCRS_sorted_filtered.bam $f)
    
#     # generate coverage metrics
#     echo "$ID generating genome coverage metrics"
#     # -a means output all genomic positions, -s means count overlapps of paried reads as one, -g means do include duplicates
#     samtools depth -a -s \
#     -g 0x400 \
#     -o "$RESULTS"/coverage/"$ID".txt \
#     "$f"

# done

echo -e "\n------------------------------------------------------------------------------\nAlignment and alignment metrics have been completed.\nIt was not triggered to end by any errors.\nOutputs can be found under:\n$RESULTS\n------------------------------------------------------------------------------\n"

## 6. VARIANT CALLING
for f in "$ALIGN"/W_A319*/*_rCRS_sorted_filtered.bam
    do
    ID=$(basename -s _rCRS_sorted_filtered.bam "$f" )

    # Run freebayes
    echo "$ID calling vars with A.Hodgkinson parameters"

    freebayes -b "$f" \
    -f "$REF_RCRS" \
    --pooled-continuous \
    --min-alternate-fraction 0.02 \
    --throw-away-complex-obs \
    -v "$VARS"/"$ID"_AH-params.vcf

    # VCF pre-processing: split VCF lines so that each line contains only one variant
    bgzip "$VARS"/"$ID"_AH-params.vcf
    tabix "$VARS"/"$ID"_AH-params.vcf.gz
    bcftools norm -m-both --output "$VARS"/"$ID"_AH-params_multiallelic.vcf "$VARS"/"$ID"_AH-params.vcf.gz

    # VCF pre-processing: force left-normalisation
    bcftools norm -f "$REF_RCRS" \
    --output "$VARS"/"$ID"_AH-params_multial-leftnorm.vcf \
    "$VARS"/"$ID"_AH-params_multiallelic.vcf

    # VCF quality control
    vcffilter -f "QUAL > 30" \
    -f "DP > 100" -f "RPL >1" -f "RPR >1" \
    "$VARS"/"$ID"_AH-params_multial-leftnorm.vcf > "$VARS"/"$ID"_AH-params_filter.vcf

    echo "$ID - converting to table"
    gatk VariantsToTable -V "$VARS"/"$ID"_AH-params_filter.vcf \
    -O  "$VARTABLE"/"$ID"_AH-params_filter_GATK_table \
    -F "POS" -F "ID" -F "REF" -F "ALT" -F "QUAL" -F "AO" -F "DP" -F "NUMALT" -F "RPL" -F "RPR" \
    --split-multi-allelic --verbosity ERROR
    
    # Clean up - remove intermediate files
    rm "$VARS"/"$ID"_AH-params_multiallelic.vcf
    rm "$VARS"/"$ID"_AH-params.vcf.gz.tbi
    rm "$VARS"/"$ID"_AH-params.vcf.gz

done