## Copyright Broad Institute, 2016
##
## This WDL pipeline implements data pre-processing and initial variant calling (GVCF
## generation) according to the GATK Best Practices (June 2016) for germline SNP and
## Indel discovery in human whole-genome sequencing (WGS) data.
##
## Requirements/expectations :
## - Human whole-genome pair-end sequencing data in unmapped BAM (uBAM) format
## - One or more read groups, one per uBAM file, all belonging to a single sample (SM)
## - Input uBAM files must additionally comply with the following requirements:
## - - filenames all have the same suffix (we use ".unmapped.bam")
## - - files must pass validation by ValidateSamFile
## - - reads are provided in query-sorted order
## - - all reads must have an RG tag
## - Reference genome must be Hg38 with ALT contigs
##
## Runtime parameters are optimized for Broad's Google Cloud Platform implementation.
## For program versions, see docker containers.
##
## LICENSING :
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.


# The internal production version of this WDL must have Private and Public block tags
# annotated to distinguish the content that can be made public versus what cannot.


# TASK DEFINITIONS


task CollectQualityYieldMetrics {
  File input_bam
  String metrics_filename
  Int disk_size
  Int preemptible_tries

  command {
    java -Xms128m -jar /usr/gitc/picard.jar \
      CollectQualityYieldMetrics \
      INPUT=${input_bam} \
      OQ=true \
      OUTPUT=${metrics_filename}
  }
  runtime {
    disks: "local-disk " + disk_size + " HDD"
    memory: "2 GB"
    preemptible: preemptible_tries
  }
  output {
    File metrics = "${metrics_filename}"
  }
}

## check the assumption that the final GVCF filename that is going to be used ends with .g.vcf.gz
task CheckFinalVcfExtension {
  String vcf_filename

  command <<<
    python <<CODE
    import os
    import sys
    filename="${vcf_filename}"
    if not filename.endswith(".vcf.gz"):
      raise Exception("input","vcf filename must end with '.g.vcf.gz', found %s"%(filename))
      sys.exit(1)
    CODE
  >>>
  runtime {
    docker: "python:2.7"
    memory: "2 GB"
  }
  output {
    String common_suffix=read_string(stdout())
  }
}


# Get version of BWA
task GetBwaVersion {
  command {
    # not setting set -o pipefail here because /bwa has a rc=1 and we dont want to allow rc=1 to succeed because
    # the sed may also fail with that error and that is something we actually want to fail on.
    /usr/gitc/bwa 2>&1 | \
    grep -e '^Version' | \
    sed 's/Version: //'
  }
  runtime {
    memory: "1 GB"
  }
  output {
    String version = read_string(stdout())
  }
}


# Read unmapped BAM, convert on-the-fly to FASTQ and stream to BWA MEM for alignment
task SamToFastqAndBwaMem {
  File input_bam
  String bwa_commandline
  String output_bam_basename
  File ref_fasta
  File ref_fasta_index
  File ref_dict

  # This is the .alt file from bwa-kit (https://github.com/lh3/bwa/tree/master/bwakit),
  # listing the reference contigs that are "alternative".
  File ref_alt

  File ref_amb
  File ref_ann
  File ref_bwt
  File ref_pac
  File ref_sa
  Int disk_size
  Int preemptible_tries

  command <<<
    set -o pipefail
    set -e

    # set the bash variable needed for the command-line
    bash_ref_fasta=${ref_fasta}
    # if ref_alt has data in it,
    if [ -s ${ref_alt} ]; then
      java -Xms3000m -jar /usr/gitc/picard.jar \
        SamToFastq \
        INPUT=${input_bam} \
        FASTQ=/dev/stdout \
        INTERLEAVE=true \
        NON_PF=true | \
      /usr/gitc/${bwa_commandline} /dev/stdin -  2> >(tee ${output_bam_basename}.bwa.stderr.log >&2) | \
      samtools view -1 - > ${output_bam_basename}.bam

      grep -m1 "read .* ALT contigs" ${output_bam_basename}.bwa.stderr.log | \
      grep -v "read 0 ALT contigs"

    # else ref_alt is empty or could not be found
    else
      exit 1;
    fi
  >>>
  runtime {
    preemptible: preemptible_tries
    memory: "14 GB"
    cpu: "16"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File output_bam = "${output_bam_basename}.bam"
    File bwa_stderr_log = "${output_bam_basename}.bwa.stderr.log"
  }
}


# Merge original input uBAM file with BWA-aligned BAM file
task MergeBamAlignment {
  File unmapped_bam
  String bwa_commandline
  String bwa_version
  File aligned_bam
  String output_bam_basename
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  Int disk_size
  Int preemptible_tries

  command {
    # set the bash variable needed for the command-line
    bash_ref_fasta=${ref_fasta}
    java -Xms3000m -jar /usr/gitc/picard.jar \
      MergeBamAlignment \
      VALIDATION_STRINGENCY=SILENT \
      EXPECTED_ORIENTATIONS=FR \
      ATTRIBUTES_TO_RETAIN=X0 \
      ALIGNED_BAM=${aligned_bam} \
      UNMAPPED_BAM=${unmapped_bam} \
      OUTPUT=${output_bam_basename}.bam \
      REFERENCE_SEQUENCE=${ref_fasta} \
      PAIRED_RUN=true \
      SORT_ORDER="unsorted" \
      IS_BISULFITE_SEQUENCE=false \
      ALIGNED_READS_ONLY=false \
      CLIP_ADAPTERS=false \
      MAX_RECORDS_IN_RAM=2000000 \
      ADD_MATE_CIGAR=true \
      MAX_INSERTIONS_OR_DELETIONS=-1 \
      PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
      PROGRAM_RECORD_ID="bwamem" \
      PROGRAM_GROUP_VERSION="${bwa_version}" \
      PROGRAM_GROUP_COMMAND_LINE="${bwa_commandline}" \
      PROGRAM_GROUP_NAME="bwamem" \
      UNMAPPED_READ_STRATEGY=COPY_TO_TAG \
      ALIGNER_PROPER_PAIR_FLAGS=true \
      UNMAP_CONTAMINANT_READS=true
  }
  runtime {
    preemptible: preemptible_tries
    memory: "3500 MB"
    cpu: "1"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File output_bam = "${output_bam_basename}.bam"
  }
}


# Sort BAM file by coordinate order and fix tag values for NM and UQ
task SortAndFixTags {
  File input_bam
  String output_bam_basename
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  Int preemptible_tries
  Int disk_size

  command {
    set -o pipefail
    java -Xms4000m -jar /usr/gitc/picard.jar \
    SortSam \
    INPUT=${input_bam} \
    OUTPUT=/dev/stdout \
    SORT_ORDER="coordinate" \
    CREATE_INDEX=false \
    CREATE_MD5_FILE=false | \
    java -Xms500m -jar /usr/gitc/picard.jar \
    SetNmAndUqTags \
    INPUT=/dev/stdin \
    OUTPUT=${output_bam_basename}.bam \
    CREATE_INDEX=true \
    CREATE_MD5_FILE=true \
    REFERENCE_SEQUENCE=${ref_fasta}
  }
  runtime {
    disks: "local-disk " + disk_size + " HDD"
    cpu: "1"
    memory: "5000 MB"
    preemptible: preemptible_tries
  }
  output {
    File output_bam = "${output_bam_basename}.bam"
    File output_bam_index = "${output_bam_basename}.bai"
    File output_bam_md5 = "${output_bam_basename}.bam.md5"
  }
}


# The tasks CollectUnsortedReadgroupBamQualityMetrics, CollectReadgroupBamQualityMetrics,
# and CollectAggregationMetrics are all variations of CollectMultipleMetrics.
# They have been separated because we cannot glob their outputs.
# When glob works again, they can be merged back into one task.


task CollectUnsortedReadgroupBamQualityMetrics {
  File input_bam
  String output_bam_prefix
  Int disk_size

  command {
    java -Xms5000m -jar /usr/gitc/picard.jar \
      CollectMultipleMetrics \
      INPUT=${input_bam} \
      OUTPUT=${output_bam_prefix} \
      ASSUME_SORTED=true \
      PROGRAM="null" \
      PROGRAM="CollectBaseDistributionByCycle" \
      PROGRAM="CollectInsertSizeMetrics" \
      PROGRAM="MeanQualityByCycle" \
      PROGRAM="QualityScoreDistribution" \
      METRIC_ACCUMULATION_LEVEL="null" \
      METRIC_ACCUMULATION_LEVEL="ALL_READS"

    touch ${output_bam_prefix}.insert_size_metrics
    touch ${output_bam_prefix}.insert_size_histogram.pdf

  }
  runtime {
    memory: "7 GB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File base_distribution_by_cycle_pdf = "${output_bam_prefix}.base_distribution_by_cycle.pdf"
    File base_distribution_by_cycle_metrics = "${output_bam_prefix}.base_distribution_by_cycle_metrics"
    File insert_size_histogram_pdf = "${output_bam_prefix}.insert_size_histogram.pdf"
    File insert_size_metrics = "${output_bam_prefix}.insert_size_metrics"
    File quality_by_cycle_pdf = "${output_bam_prefix}.quality_by_cycle.pdf"
    File quality_by_cycle_metrics = "${output_bam_prefix}.quality_by_cycle_metrics"
    File quality_distribution_pdf = "${output_bam_prefix}.quality_distribution.pdf"
    File quality_distribution_metrics = "${output_bam_prefix}.quality_distribution_metrics"
  }
}


task CollectReadgroupBamQualityMetrics {
  File input_bam
  File input_bam_index
  String output_bam_prefix
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  Int disk_size

  command {
    java -Xms5000m -jar /usr/gitc/picard.jar \
      CollectMultipleMetrics \
      INPUT=${input_bam} \
      REFERENCE_SEQUENCE=${ref_fasta} \
      OUTPUT=${output_bam_prefix} \
      ASSUME_SORTED=true \
      PROGRAM="null" \
      PROGRAM="CollectAlignmentSummaryMetrics" \
      PROGRAM="CollectGcBiasMetrics" \
      METRIC_ACCUMULATION_LEVEL="null" \
      METRIC_ACCUMULATION_LEVEL="READ_GROUP"
  }
  runtime {
    memory: "7 GB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File alignment_summary_metrics = "${output_bam_prefix}.alignment_summary_metrics"
    File gc_bias_detail_metrics = "${output_bam_prefix}.gc_bias.detail_metrics"
    File gc_bias_pdf = "${output_bam_prefix}.gc_bias.pdf"
    File gc_bias_summary_metrics = "${output_bam_prefix}.gc_bias.summary_metrics"
  }
}


task CollectAggregationMetrics {
  File input_bam
  File input_bam_index
  String output_bam_prefix
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  Int disk_size

  command {
    java -Xms5000m -jar /usr/gitc/picard.jar \
      CollectMultipleMetrics \
      INPUT=${input_bam} \
      REFERENCE_SEQUENCE=${ref_fasta} \
      OUTPUT=${output_bam_prefix} \
      ASSUME_SORTED=true \
      PROGRAM="null" \
      PROGRAM="CollectAlignmentSummaryMetrics" \
      PROGRAM="CollectInsertSizeMetrics" \
      PROGRAM="CollectSequencingArtifactMetrics" \
      PROGRAM="CollectGcBiasMetrics" \
      PROGRAM="QualityScoreDistribution" \
      METRIC_ACCUMULATION_LEVEL="null" \
      METRIC_ACCUMULATION_LEVEL="SAMPLE" \
      METRIC_ACCUMULATION_LEVEL="LIBRARY"

    touch ${output_bam_prefix}.insert_size_metrics
    touch ${output_bam_prefix}.insert_size_histogram.pdf
  }
  runtime {
    memory: "7 GB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File alignment_summary_metrics = "${output_bam_prefix}.alignment_summary_metrics"
    File bait_bias_detail_metrics = "${output_bam_prefix}.bait_bias_detail_metrics"
    File bait_bias_summary_metrics = "${output_bam_prefix}.bait_bias_summary_metrics"
    File gc_bias_detail_metrics = "${output_bam_prefix}.gc_bias.detail_metrics"
    File gc_bias_pdf = "${output_bam_prefix}.gc_bias.pdf"
    File gc_bias_summary_metrics = "${output_bam_prefix}.gc_bias.summary_metrics"
    File insert_size_histogram_pdf = "${output_bam_prefix}.insert_size_histogram.pdf"
    File insert_size_metrics = "${output_bam_prefix}.insert_size_metrics"
    File pre_adapter_detail_metrics = "${output_bam_prefix}.pre_adapter_detail_metrics"
    File pre_adapter_summary_metrics = "${output_bam_prefix}.pre_adapter_summary_metrics"
    File quality_distribution_pdf = "${output_bam_prefix}.quality_distribution.pdf"
    File quality_distribution_metrics = "${output_bam_prefix}.quality_distribution_metrics"
  }
}


# Mark duplicate reads to avoid counting non-independent observations
task MarkDuplicates {
  Array[File] input_bams
  String output_bam_basename
  String metrics_filename
  Int disk_size
  # The program default for READ_NAME_REGEX is appropriate in nearly every case.
  # Sometimes we wish to supply "null" in order to turn off optical duplicate detection
  # This can be desirable if you don't mind the estimated library size being wrong and optical duplicate detection is taking >7 days and failing
  String? read_name_regex

 # Task is assuming query-sorted input so that the Secondary and Supplementary reads get marked correctly
 # This works because the output of BWA is query-grouped, and thus so is the output of MergeBamAlignment.
 # While query-grouped isn't actually query-sorted, it's good enough for MarkDuplicates with ASSUME_SORT_ORDER="queryname"
  command {
    java -Xms4000m -jar /usr/gitc/picard.jar \
      MarkDuplicates \
      INPUT=${sep=' INPUT=' input_bams} \
      OUTPUT=${output_bam_basename}.bam \
      METRICS_FILE=${metrics_filename} \
      VALIDATION_STRINGENCY=SILENT \
      ${"READ_NAME_REGEX=" + read_name_regex} \
      OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
      ASSUME_SORT_ORDER="queryname"
      CREATE_MD5_FILE=true
  }
  runtime {
    memory: "7 GB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File output_bam = "${output_bam_basename}.bam"
    File duplicate_metrics = "${metrics_filename}"
  }
}


# Generate sets of intervals for scatter-gathering over chromosomes
task CreateSequenceGroupingTSV {
  File ref_dict
  Int preemptible_tries

  # Use python to create the Sequencing Groupings used for BQSR and PrintReads Scatter.  It outputs to stdout
  # where it is parsed into a wdl Array[Array[String]]
  # e.g. [["1"], ["2"], ["3", "4"], ["5"], ["6", "7", "8"]]
  command <<<
    python <<CODE
    with open("${ref_dict}", "r") as ref_dict_file:
        sequence_tuple_list = []
        longest_sequence = 0
        for line in ref_dict_file:
            if line.startswith("@SQ"):
                line_split = line.split("\t")
                # (Sequence_Name, Sequence_Length)
                sequence_tuple_list.append((line_split[1].split("SN:")[1], int(line_split[2].split("LN:")[1])))
        longest_sequence = sorted(sequence_tuple_list, key=lambda x: x[1], reverse=True)[0][1]
    # We are adding this to the intervals because hg38 has contigs named with embedded colons and a bug in GATK strips off
    # the last element after a :, so we add this as a sacrificial element.
    hg38_protection_tag = ":1+"
    # initialize the tsv string with the first sequence
    tsv_string = sequence_tuple_list[0][0] + hg38_protection_tag
    temp_size = sequence_tuple_list[0][1]
    for sequence_tuple in sequence_tuple_list[1:]:
        if temp_size + sequence_tuple[1] <= longest_sequence:
            temp_size += sequence_tuple[1]
            tsv_string += "\t" + sequence_tuple[0] + hg38_protection_tag
        else:
            tsv_string += "\n" + sequence_tuple[0] + hg38_protection_tag
            temp_size = sequence_tuple[1]
    # add the unmapped sequences as a separate line to ensure that they are recalibrated as well
    with open("sequence_grouping.txt","w") as tsv_file:
      tsv_file.write(tsv_string)
      tsv_file.close()

    tsv_string += '\n' + "unmapped"

    with open("sequence_grouping_with_unmapped.txt","w") as tsv_file_with_unmapped:
      tsv_file_with_unmapped.write(tsv_string)
      tsv_file_with_unmapped.close()
    CODE
  >>>
  runtime {
    preemptible: preemptible_tries
    docker: "python:2.7"
    memory: "2 GB"
  }
  output {
    Array[Array[String]] sequence_grouping = read_tsv("sequence_grouping.txt")
    Array[Array[String]] sequence_grouping_with_unmapped = read_tsv("sequence_grouping_with_unmapped.txt")
  }
}


# Generate Base Quality Score Recalibration (BQSR) model
task BaseRecalibrator {
  File input_bam
  File input_bam_index
  String recalibration_report_filename
  Array[String] sequence_group_interval
  File dbSNP_vcf
  File dbSNP_vcf_index
  Array[File] known_indels_sites_VCFs
  Array[File] known_indels_sites_indices
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  Int disk_size
  Int preemptible_tries

  command {
    java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal \
      -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails \
      -Xloggc:gc_log.log -Dsamjdk.use_async_io=false -Xms4000m \
      -jar /usr/gitc/GATK4.jar \
      BaseRecalibrator \
      -R ${ref_fasta} \
      -I ${input_bam} \
      --useOriginalQualities \
      -O ${recalibration_report_filename} \
      -knownSites ${dbSNP_vcf} \
      -knownSites ${sep=" -knownSites " known_indels_sites_VCFs} \
      -L ${sep=" -L " sequence_group_interval}
  }
  runtime {
    preemptible: preemptible_tries
    memory: "6 GB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File recalibration_report = "${recalibration_report_filename}"
    #this output is only for GOTC STAGING to give some GC statistics to the GATK4 team
    #File gc_logs = "gc_log.log"
  }
}


# Apply Base Quality Score Recalibration (BQSR) model
task ApplyBQSR {
  File input_bam
  File input_bam_index
  String output_bam_basename
  File recalibration_report
  Array[String] sequence_group_interval
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  Int disk_size
  Int preemptible_tries

  command {
    java -XX:+PrintFlagsFinal -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps \
      -XX:+PrintGCDetails -Xloggc:gc_log.log -Dsamjdk.use_async_io=false \
      -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms3000m \
      -jar /usr/gitc/GATK4.jar \
      ApplyBQSR \
      --createOutputBamMD5 \
      --addOutputSAMProgramRecord \
      -R ${ref_fasta} \
      -I ${input_bam} \
      --useOriginalQualities \
      -O ${output_bam_basename}.bam \
      -bqsr ${recalibration_report} \
      -SQQ 10 -SQQ 20 -SQQ 30 \
      -L ${sep=" -L " sequence_group_interval}
  }
  runtime {
    preemptible: preemptible_tries
    memory: "3500 MB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File recalibrated_bam = "${output_bam_basename}.bam"
    File recalibrated_bam_checksum = "${output_bam_basename}.bam.md5"
    #this output is only for GOTC STAGING to give some GC statistics to the GATK4 team
    #File gc_logs = "gc_log.log"
  }
}


# Combine multiple recalibration tables from scattered BaseRecalibrator runs
task GatherBqsrReports {
  Array[File] input_bqsr_reports
  String output_report_filename
  Int disk_size
  Int preemptible_tries

  command {
    java -Xms3000m -jar /usr/gitc/GATK4.jar \
      GatherBQSRReports \
      -I ${sep=' -I ' input_bqsr_reports} \
      -O ${output_report_filename}
    }
  runtime {
    preemptible: preemptible_tries
    memory: "3500 MB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File output_bqsr_report = "${output_report_filename}"
  }
}


# Combine multiple recalibrated BAM files from scattered ApplyRecalibration runs
task GatherBamFiles {
  Array[File] input_bams
  String output_bam_basename
  Int disk_size
  Int preemptible_tries

  command {
    java -Xms2000m -jar /usr/gitc/picard.jar \
      GatherBamFiles \
      INPUT=${sep=' INPUT=' input_bams} \
      OUTPUT=${output_bam_basename}.bam \
      CREATE_INDEX=true \
      CREATE_MD5_FILE=true
    }
  runtime {
    preemptible: preemptible_tries
    memory: "3 GB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File output_bam = "${output_bam_basename}.bam"
    File output_bam_index = "${output_bam_basename}.bai"
    File output_bam_md5 = "${output_bam_basename}.bam.md5"
  }
}


task ValidateSamFile {
  File input_bam
  File? input_bam_index
  String report_filename
  File? ref_dict
  File? ref_fasta
  File? ref_fasta_index
  Int? max_output
  Array[String]? ignore
  Int disk_size
  Int preemptible_tries

  command {
    java -Xms4000m -jar /usr/gitc/picard.jar \
      ValidateSamFile \
      INPUT=${input_bam} \
      OUTPUT=${report_filename} \
      ${"REFERENCE_SEQUENCE=" + ref_fasta} \
      ${"MAX_OUTPUT=" + max_output} \
      IGNORE=${default="null" sep=" IGNORE=" ignore} \
      MODE=VERBOSE \
      IS_BISULFITE_SEQUENCED=false
  }
  runtime {
    preemptible: preemptible_tries
    memory: "7 GB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File report = "${report_filename}"
  }
}


# TODO - CollectWgsMetrics and CollectRawWgsMetrics can be combined by a locus-iterating collector
task CollectWgsMetrics {
  File input_bam
  File input_bam_index
  String metrics_filename
  File wgs_coverage_interval_list
  File ref_fasta
  File ref_fasta_index
  Int disk_size

  command {
    java -Xms2000m -jar /usr/gitc/picard.jar \
      CollectWgsMetrics \
      INPUT=${input_bam} \
      VALIDATION_STRINGENCY=SILENT \
      REFERENCE_SEQUENCE=${ref_fasta} \
      INTERVALS=${wgs_coverage_interval_list} \
      OUTPUT=${metrics_filename}
  }
  runtime {
    memory: "3 GB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File metrics = "${metrics_filename}"
  }
}


task CollectRawWgsMetrics {
  File input_bam
  File input_bam_index
  String metrics_filename
  File wgs_coverage_interval_list
  File ref_fasta
  File ref_fasta_index
  Int disk_size

  command {
    java -Xms2000m -jar /usr/gitc/picard.jar \
      CollectRawWgsMetrics \
      INPUT=${input_bam} \
      VALIDATION_STRINGENCY=SILENT \
      REFERENCE_SEQUENCE=${ref_fasta} \
      INTERVALS=${wgs_coverage_interval_list} \
      OUTPUT=${metrics_filename}
  }
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:2.2.5-1486412288"
    memory: "3 GB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File metrics = "${metrics_filename}"
  }
}


task CalculateReadGroupChecksum {
  File input_bam
  File input_bam_index
  String read_group_md5_filename
  Int disk_size
  Int preemptible_tries

  command {
    java -Xms1000m -jar /usr/gitc/picard.jar \
      CalculateReadGroupChecksum \
      INPUT=${input_bam} \
      OUTPUT=${read_group_md5_filename}
  }
  runtime {
    preemptible: preemptible_tries
    memory: "2 GB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File md5_file = "${read_group_md5_filename}"
  }
}


# Notes on the contamination estimate:
# The contamination value is read from the FREEMIX field of the selfSM file output by verifyBamId
#
# In Zamboni production, this value is stored directly in METRICS.AGGREGATION_CONTAM
#
# Contamination is also stored in GVCF_CALLING and thereby passed to HAPLOTYPE_CALLER
# But first, it is divided by an underestimation factor thusly:
#   float(FREEMIX) / ContaminationUnderestimationFactor
#     where the denominator is hardcoded in Zamboni:
#     val ContaminationUnderestimationFactor = 0.75f
#
# Here, I am handling this by returning both the original selfSM file for reporting, and the adjusted
#   contamination estimate for use in variant calling
task CheckContamination {
  File input_bam
  File input_bam_index
  File contamination_sites_vcf
  File contamination_sites_vcf_index
  String output_prefix
  Int disk_size
  Int preemptible_tries

  # Having to do this as a 2-step command in heredoc syntax, adding a python step to read the metrics
  # This is a hack until read_object() is supported by Cromwell.
  # It relies on knowing that there is only one data row in the 2-row selfSM TSV file
  # Piping output of verifyBamId to /dev/null so only stdout is from the python command
  command <<<
    set -e

    /usr/gitc/verifyBamID \
    --verbose \
    --ignoreRG \
    --vcf ${contamination_sites_vcf} \
    --out ${output_prefix} \
    --bam ${input_bam} \
    1>/dev/null

    python3 <<CODE
    import csv
    import sys
    with open('${output_prefix}.selfSM') as selfSM:
      reader = csv.DictReader(selfSM, delimiter='\t')
      i = 0
      for row in reader:
        if float(row["FREELK0"])==0 and float(row["FREELK1"])==0:
    # a zero value for the likelihoods implies no data. This usually indicates a problem rather than a real event.
    # if the bam isn't really empty, this is probably due to the use of a incompatible reference build between
    # vcf and bam.
          sys.stderr.write("Found zero likelihoods. Bam is either very-very shallow, or aligned to the wrong reference (relative to the vcf).")
          sys.exit(1)
        print(float(row["FREEMIX"])/0.75)
        i = i + 1
    # there should be exactly one row, and if this isn't the case the format of the output is unexpectedly different
    # and the results are not reliable.
        if i != 1:
          sys.stderr.write("Found %d rows in .selfSM file. Was expecting exactly 1. This is an error"%(i))
          sys.exit(2)
    CODE
  >>>
  runtime {
    preemptible: preemptible_tries
    memory: "2 GB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File selfSM = "${output_prefix}.selfSM"
    File depthSM = "${output_prefix}.depthSM"
    File log = "${output_prefix}.log"

    # I would like to do the following, however:
    # The object is read as a string
    # explicit string->float coercion via float(), as shown below, is supported by Cromwell
    # the interim value cannot be stored as a string and then assigned to a float. Variables intialized in output cannot be dereferenced in output.
    # Float contamination = float(read_object(${output_prefix} + ".selfSM").FREEMIX) / 0.75

    # In the interim, get the value from the python hack above:
    Float contamination = read_float(stdout())
  }
}

# This task calls picard's IntervalListTools to scatter the input interval list into scatter_count sub interval lists
# Note that the number of sub interval lists may not be exactly equal to scatter_count.  There may be slightly more or less.
# Thus we have the block of python to count the number of generated sub interval lists.
task ScatterIntervalList {
  File interval_list
  Int scatter_count
  Int break_bands_at_multiples_of

  command <<<
    set -e
    mkdir out
    java -Xms1g -jar /usr/gitc/picard.jar \
    IntervalListTools \
    SCATTER_COUNT=${scatter_count} \
    SUBDIVISION_MODE=BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \
    UNIQUE=true \
    SORT=true \
    BREAK_BANDS_AT_MULTIPLES_OF=${break_bands_at_multiples_of} \
    INPUT=${interval_list} \
    OUTPUT=out
    python3 <<CODE
    import glob, os
    # Works around a JES limitation where multiples files with the same name overwrite each other when globbed
    intervals = sorted(glob.glob("out/*/*.interval_list"))
    for i, interval in enumerate(intervals):
      (directory, filename) = os.path.split(interval)
      newName = os.path.join(directory, str(i + 1) + filename)
      os.rename(interval, newName)
    print(len(intervals))
    CODE
  >>>
  output {
    Array[File] out = glob("out/*/*.interval_list")
    Int interval_count = read_int(stdout())
  }

  runtime {
    memory: "2 GB"
  }
}



# Call variants on a single sample with HaplotypeCaller to produce a GVCF
task HaplotypeCaller {
  File input_bam
  File input_bam_index
  File interval_list
  String gvcf_basename
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  Float? contamination
  Int disk_size
  Int preemptible_tries

  # tried to find lowest memory variable where it would still work, might change once tested on JES
  command {
    java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms8000m \
      -jar /usr/gitc/GATK35.jar \
      -T HaplotypeCaller \
      -R ${ref_fasta} \
      -o ${gvcf_basename}.vcf.gz \
      -I ${input_bam} \
      -L ${interval_list} \
      -ERC GVCF \
      --max_alternate_alleles 3 \
      -variant_index_parameter 128000 \
      -variant_index_type LINEAR \
      -contamination ${default=0 contamination} \
      --read_filter OverclippedRead
  }
  runtime {
    preemptible: preemptible_tries
    memory: "10 GB"
    cpu: "1"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File output_gvcf = "${gvcf_basename}.vcf.gz"
    File output_gvcf_index = "${gvcf_basename}.vcf.gz.tbi"
  }
}


# Combine multiple VCFs or GVCFs from scattered HaplotypeCaller runs
task MergeVCFs {
  Array[File] input_vcfs
  Array[File] input_vcfs_indexes
  String output_vcf_name
  Int disk_size
  Int preemptible_tries

  # using MergeVcfs instead of GatherVcfs so we can create indices
  # WARNING 2015-10-28 15:01:48 GatherVcfs  Index creation not currently supported when gathering block compressed VCFs.
  command {
    java -Xms2000m -jar /usr/gitc/picard.jar \
    MergeVcfs \
    INPUT=${sep=' INPUT=' input_vcfs} \
    OUTPUT=${output_vcf_name}
  }
  output {
    File output_vcf = "${output_vcf_name}"
    File output_vcf_index = "${output_vcf_name}.tbi"
  }
  runtime {
    preemptible: preemptible_tries
    memory: "3 GB"
    disks: "local-disk " + disk_size + " HDD"
  }
}


task ValidateGVCF {
  File input_vcf
  File input_vcf_index
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  File dbSNP_vcf
  File dbSNP_vcf_index
  File wgs_calling_interval_list
  Int disk_size
  Int preemptible_tries

  command {
    java -Xms8000m -jar /usr/gitc/GATK36.jar \
    -T ValidateVariants \
    -V ${input_vcf} \
    -R ${ref_fasta} \
    -L ${wgs_calling_interval_list} \
    -gvcf \
    --validationTypeToExclude ALLELES \
    --reference_window_stop 208 -U  \
    --dbsnp ${dbSNP_vcf}
  }
  runtime {
    preemptible: preemptible_tries
    memory: "10 GB"
    disks: "local-disk " + disk_size + " HDD"
  }
}


task CollectGvcfCallingMetrics {
  File input_vcf
  File input_vcf_index
  String metrics_basename
  File dbSNP_vcf
  File dbSNP_vcf_index
  File ref_dict
  File wgs_evaluation_interval_list
  Int disk_size
  Int preemptible_tries

  command {
    java -Xms2000m -jar /usr/gitc/picard.jar \
      CollectVariantCallingMetrics \
      INPUT=${input_vcf} \
      OUTPUT=${metrics_basename} \
      DBSNP=${dbSNP_vcf} \
      SEQUENCE_DICTIONARY=${ref_dict} \
      TARGET_INTERVALS=${wgs_evaluation_interval_list} \
      GVCF_INPUT=true
  }
  runtime {
    preemptible: preemptible_tries
    memory: "3 GB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File summary_metrics = "${metrics_basename}.variant_calling_summary_metrics"
    File detail_metrics = "${metrics_basename}.variant_calling_detail_metrics"
  }
}


# Convert BAM file to CRAM format
task ConvertToCram {
  File input_bam
  File ref_fasta
  File ref_fasta_index
  String output_basename
  Int disk_size
  Int preemptible_tries

  command <<<
    set -e
    set -o pipefail

    samtools view -C -T ${ref_fasta} ${input_bam} | \
    tee ${output_basename}.cram | \
    md5sum | awk '{print $1}' > ${output_basename}.cram.md5

    # Create REF_CACHE. Used when indexing a CRAM
    seq_cache_populate.pl -root ./ref/cache ${ref_fasta}
    export REF_PATH=:
    export REF_CACHE=./ref/cache/%2s/%2s/%s

    samtools index ${output_basename}.cram
    mv ${output_basename}.cram.crai ${output_basename}.crai
  >>>
  runtime {
    preemptible: preemptible_tries
    memory: "3 GB"
    cpu: "1"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File output_cram = "${output_basename}.cram"
    File output_cram_index = "${output_basename}.crai"
    File output_cram_md5 = "${output_basename}.cram.md5"
  }
}

task CramToBam {
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  File cram_file
  String output_basename

  Int disk_size

command <<<
  set -e
  set -o pipefail

  samtools view -h -T ${ref_fasta} ${cram_file} |
  samtools view -b -o ${output_basename}.bam -
  samtools index -b ${output_basename}.bam
  mv ${output_basename}.bam.bai ${output_basename}.bai
  >>>
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:2.2.5-1486412288"
    memory: "3 GB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File output_bam = "${output_basename}.bam"
    File output_bam_index = "${output_basename}.bai"
  }
}


# WORKFLOW DEFINITION
workflow PairedEndSingleSampleWorkflow {

  
  File contamination_sites_vcf
  File contamination_sites_vcf_index
  File wgs_evaluation_interval_list
  File wgs_coverage_interval_list

  
  String sample_name
  String base_file_name
  String final_gvcf_name
  Array[File] flowcell_unmapped_bams
  String unmapped_bam_suffix

  File wgs_calling_interval_list
  Int haplotype_scatter_count
  Int break_bands_at_multiples_of

  File ref_fasta
  File ref_fasta_index
  File ref_dict
  File ref_alt
  File ref_bwt
  File ref_sa
  File ref_amb
  File ref_ann
  File ref_pac

  File dbSNP_vcf
  File dbSNP_vcf_index
  Array[File] known_indels_sites_VCFs
  Array[File] known_indels_sites_indices

  Int flowcell_small_disk
  Int flowcell_medium_disk
  Int agg_small_disk
  Int agg_medium_disk
  Int agg_large_disk
  Int preemptible_tries
  Int agg_preemptible_tries

  String bwa_commandline="bwa mem -K 100000000 -p -v 3 -t 16 -Y $bash_ref_fasta"

  String recalibrated_bam_basename = base_file_name + ".aligned.duplicates_marked.recalibrated"

  
  # Get the version of BWA to include in the PG record in the header of the BAM produced
  # by MergeBamAlignment.
  call GetBwaVersion

  call CheckFinalVcfExtension {
     input:
        vcf_filename = final_gvcf_name
   }

  
  # Align flowcell-level unmapped input bams in parallel
  scatter (unmapped_bam in flowcell_unmapped_bams) {

    String sub_strip_path = "gs://.*/"
    String sub_strip_unmapped = unmapped_bam_suffix + "$"
    String sub_sub = sub(sub(unmapped_bam, sub_strip_path, ""), sub_strip_unmapped, "")
    
    call CollectQualityYieldMetrics {
      input:
        input_bam = unmapped_bam,
        metrics_filename = sub_sub + ".unmapped.quality_yield_metrics",
        disk_size = flowcell_small_disk,
        preemptible_tries = preemptible_tries
    }

    
    # Map reads to reference
    call SamToFastqAndBwaMem {
      input:
        input_bam = unmapped_bam,
        bwa_commandline = bwa_commandline,
        output_bam_basename = sub_sub + ".unmerged",
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        ref_alt = ref_alt,
        ref_bwt = ref_bwt,
        ref_amb = ref_amb,
        ref_ann = ref_ann,
        ref_pac = ref_pac,
        ref_sa = ref_sa,
        disk_size = flowcell_medium_disk,
        preemptible_tries = preemptible_tries
     }

    
    # Merge original uBAM and BWA-aligned BAM
    call MergeBamAlignment {
      input:
        unmapped_bam = unmapped_bam,
        bwa_commandline = bwa_commandline,
        bwa_version = GetBwaVersion.version,
        aligned_bam = SamToFastqAndBwaMem.output_bam,
        output_bam_basename = sub_sub + ".aligned.unsorted",
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        disk_size = flowcell_medium_disk,
        preemptible_tries = preemptible_tries
    }

    
    # Sort and fix tags in the merged BAM
    call SortAndFixTags as SortAndFixReadGroupBam {
      input:
      input_bam = MergeBamAlignment.output_bam,
      output_bam_basename = sub_sub + ".sorted",
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      preemptible_tries = preemptible_tries,
      disk_size = flowcell_medium_disk
    }

    
    # called to help in finding problems early.
    # if too time consuming and not helpful, can be removed.
    call ValidateSamFile as ValidateReadGroupSamFile {
      input:
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        input_bam = SortAndFixReadGroupBam.output_bam,
        report_filename = sub_sub + ".validation_report",
        disk_size = flowcell_medium_disk,
        preemptible_tries = preemptible_tries
    }

    
    # no reference as the input here is unsorted, providing a reference would cause an error
    call CollectUnsortedReadgroupBamQualityMetrics {
      input:
        input_bam = MergeBamAlignment.output_bam,
        output_bam_prefix = sub_sub + ".readgroup",
        disk_size = flowcell_medium_disk
    }
  
  }

  
  # Aggregate aligned+merged flowcell BAM files and mark duplicates
  call MarkDuplicates {
    input:
      input_bams = MergeBamAlignment.output_bam,
      output_bam_basename = base_file_name + ".aligned.unsorted.duplicates_marked",
      metrics_filename = base_file_name + ".duplicate_metrics",
      disk_size = agg_large_disk
  }

  
  # Sort aggregated+deduped BAM file and fix tags
  call SortAndFixTags as SortAndFixSampleBam {
    input:
      input_bam = MarkDuplicates.output_bam,
      output_bam_basename = base_file_name + ".aligned.duplicate_marked.sorted",
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      disk_size = agg_large_disk,
      preemptible_tries = 0
  }


  
  # Create list of sequences for scatter-gather parallelization
  call CreateSequenceGroupingTSV {
    input:
      ref_dict = ref_dict,
      preemptible_tries = preemptible_tries
  }

  
  call CheckContamination {
    input:
      input_bam = SortAndFixSampleBam.output_bam,
      input_bam_index = SortAndFixSampleBam.output_bam_index,
      contamination_sites_vcf = contamination_sites_vcf,
      contamination_sites_vcf_index = contamination_sites_vcf_index,
      output_prefix = base_file_name + ".preBqsr",
      disk_size = agg_small_disk,
      preemptible_tries = agg_preemptible_tries
  }

  
  # Perform Base Quality Score Recalibration (BQSR) on the sorted BAM in parallel
  scatter (subgroup in CreateSequenceGroupingTSV.sequence_grouping) {
    # Generate the recalibration model by interval
    call BaseRecalibrator {
      input:
        input_bam = SortAndFixSampleBam.output_bam,
        input_bam_index = SortAndFixSampleBam.output_bam_index,
        recalibration_report_filename = base_file_name + ".recal_data.csv",
        sequence_group_interval = subgroup,
        dbSNP_vcf = dbSNP_vcf,
        dbSNP_vcf_index = dbSNP_vcf_index,
        known_indels_sites_VCFs = known_indels_sites_VCFs,
        known_indels_sites_indices = known_indels_sites_indices,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        disk_size = agg_small_disk,
        preemptible_tries = agg_preemptible_tries
    }
  }

  
  # Merge the recalibration reports resulting from by-interval recalibration
  call GatherBqsrReports {
    input:
      input_bqsr_reports = BaseRecalibrator.recalibration_report,
      output_report_filename = base_file_name + ".recal_data.csv",
      disk_size = flowcell_small_disk,
      preemptible_tries = preemptible_tries
  }

  scatter (subgroup in CreateSequenceGroupingTSV.sequence_grouping_with_unmapped) {
  
    # Apply the recalibration model by interval
    call ApplyBQSR {
      input:
        input_bam = SortAndFixSampleBam.output_bam,
        input_bam_index = SortAndFixSampleBam.output_bam_index,
        output_bam_basename = recalibrated_bam_basename,
        recalibration_report = GatherBqsrReports.output_bqsr_report,
        sequence_group_interval = subgroup,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        disk_size = agg_small_disk,
        preemptible_tries = agg_preemptible_tries
    }
  }
  
  # Merge the recalibrated BAM files resulting from by-interval recalibration

  call GatherBamFiles {
    input:
      input_bams = ApplyBQSR.recalibrated_bam,
      output_bam_basename = base_file_name,
      disk_size = agg_large_disk,
      preemptible_tries = agg_preemptible_tries
  }

  
  call CollectReadgroupBamQualityMetrics {
    input:
      input_bam = GatherBamFiles.output_bam,
      input_bam_index = GatherBamFiles.output_bam_index,
      output_bam_prefix = base_file_name + ".readgroup",
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      disk_size = agg_small_disk
  }

  
  call ValidateSamFile as ValidateAggregatedSamFile {
    input:
      input_bam = GatherBamFiles.output_bam,
      input_bam_index = GatherBamFiles.output_bam_index,
      report_filename = base_file_name + ".validation_report",
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      disk_size = agg_small_disk,
      preemptible_tries = agg_preemptible_tries
  }

  
  call CollectAggregationMetrics {
    input:
      input_bam = GatherBamFiles.output_bam,
      input_bam_index = GatherBamFiles.output_bam_index,
      output_bam_prefix = base_file_name,
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      disk_size = agg_small_disk
  }

  
  call CollectWgsMetrics {
    input:
      input_bam = GatherBamFiles.output_bam,
      input_bam_index = GatherBamFiles.output_bam_index,
      metrics_filename = base_file_name + ".wgs_metrics",
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      wgs_coverage_interval_list = wgs_coverage_interval_list,
      disk_size = agg_small_disk
  }

  
  call CollectRawWgsMetrics {
    input:
      input_bam = GatherBamFiles.output_bam,
      input_bam_index = GatherBamFiles.output_bam_index,
      metrics_filename = base_file_name + ".raw_wgs_metrics",
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      wgs_coverage_interval_list = wgs_coverage_interval_list,
      disk_size = agg_small_disk
  }

  
  call CalculateReadGroupChecksum {
    input:
      input_bam = GatherBamFiles.output_bam,
      input_bam_index = GatherBamFiles.output_bam_index,
      read_group_md5_filename = recalibrated_bam_basename + ".bam.read_group_md5",
      disk_size = agg_small_disk,
      preemptible_tries = agg_preemptible_tries
  }

  
  # Convert the final merged recalibrated BAM file to CRAM format
  call ConvertToCram {
    input:
      input_bam = GatherBamFiles.output_bam,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      output_basename = base_file_name,
      disk_size = agg_medium_disk,
      preemptible_tries = agg_preemptible_tries
  }

  
# call ValidateSamFile as ValidateCram {
#   input:
#     input_bam = ConvertToCram.output_cram,
#     input_bam_index = ConvertToCram.output_cram_index,
#     report_filename = base_file_name + ".cram.validation_report",
#     ref_dict = ref_dict,
#     ref_fasta = ref_fasta,
#     ref_fasta_index = ref_fasta_index,
#     ignore = ["MISSING_TAG_NM"],
#     max_output = 1000000000,
#     disk_size = agg_small_disk,
#     preemptible_tries = agg_preemptible_tries
# }


   

  call CramToBam {
    input:
      ref_fasta = ref_fasta,
      ref_dict = ref_dict,
      ref_fasta_index = ref_fasta_index,
      cram_file = ConvertToCram.output_cram,
      output_basename = base_file_name + ".roundtrip",
      disk_size = agg_medium_disk

  }

  call ValidateSamFile as ValidateBamFromCram {
    input:
      input_bam = CramToBam.output_bam,
      input_bam_index = CramToBam.output_bam_index,
      report_filename = base_file_name + ".bam.roundtrip.validation_report",
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      max_output = 1000000000,
      disk_size = agg_small_disk,
      ignore = ["null"],
      preemptible_tries = agg_preemptible_tries
  }

  # Break the calling interval_list into sub-intervals
  # Perform variant calling on the sub-intervals, and then gather the results
  call ScatterIntervalList {
    input:
      interval_list = wgs_calling_interval_list,
      scatter_count = haplotype_scatter_count,
      break_bands_at_multiples_of = break_bands_at_multiples_of
  }

  
  # Call variants in parallel over WGS calling intervals
  scatter (index in range(ScatterIntervalList.interval_count)) {
    
    # Generate GVCF by interval
    call HaplotypeCaller {
      input:
        
        contamination = CheckContamination.contamination,
        
        input_bam = GatherBamFiles.output_bam,
        input_bam_index = GatherBamFiles.output_bam_index,
        interval_list = ScatterIntervalList.out[index],
        gvcf_basename = sample_name,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        disk_size = agg_small_disk,
        preemptible_tries = agg_preemptible_tries
     }
  }

  
  # Combine by-interval GVCFs into a single sample GVCF file
  call MergeVCFs {
    input:
      input_vcfs = HaplotypeCaller.output_gvcf,
      input_vcfs_indexes = HaplotypeCaller.output_gvcf_index,
      output_vcf_name = final_gvcf_name,
      disk_size = agg_small_disk,
      preemptible_tries = agg_preemptible_tries
  }

  
  call ValidateGVCF {
    input:
      input_vcf = MergeVCFs.output_vcf,
      input_vcf_index = MergeVCFs.output_vcf_index,
      dbSNP_vcf = dbSNP_vcf,
      dbSNP_vcf_index = dbSNP_vcf_index,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      wgs_calling_interval_list = wgs_calling_interval_list,
      disk_size = agg_small_disk,
      preemptible_tries = agg_preemptible_tries
  }

  
  call CollectGvcfCallingMetrics {
    input:
      input_vcf = MergeVCFs.output_vcf,
      input_vcf_index = MergeVCFs.output_vcf_index,
      metrics_basename = base_file_name,
      dbSNP_vcf = dbSNP_vcf,
      dbSNP_vcf_index = dbSNP_vcf_index,
      ref_dict = ref_dict,
      wgs_evaluation_interval_list = wgs_evaluation_interval_list,
      disk_size = agg_small_disk,
      preemptible_tries = agg_preemptible_tries
  }

  
  # NOTE TO COMMS: Outputs must be redacted manually pending bug fix!!!
  # Outputs that will be retained when execution is complete

  output {
    Array[File] quality_yield_metrics = CollectQualityYieldMetrics.metrics
    Array[File] validate_group_read_sam_file_out = ValidateReadGroupSamFile.report

    Array[File] unsorted_read_group_base_distribution_by_cycle_pdf = CollectUnsortedReadgroupBamQualityMetrics.base_distribution_by_cycle_pdf
    Array[File] unsorted_read_group_base_distribution_by_cycle_metrics = CollectUnsortedReadgroupBamQualityMetrics.base_distribution_by_cycle_metrics
    Array[File] unsorted_read_group_insert_size_histogram_pdf = CollectUnsortedReadgroupBamQualityMetrics.insert_size_histogram_pdf
    Array[File] unsorted_read_group_insert_size_metrics = CollectUnsortedReadgroupBamQualityMetrics.insert_size_metrics
    Array[File] unsorted_read_group_quality_by_cycle_pdf = CollectUnsortedReadgroupBamQualityMetrics.quality_by_cycle_pdf
    Array[File] unsorted_read_group_quality_by_cycle_metrics = CollectUnsortedReadgroupBamQualityMetrics.quality_by_cycle_metrics
    Array[File] unsorted_read_group_quality_distribution_pdf = CollectUnsortedReadgroupBamQualityMetrics.quality_distribution_pdf
    Array[File] unsorted_read_group_quality_distribution_metrics = CollectUnsortedReadgroupBamQualityMetrics.quality_distribution_metrics

    File read_group_alignment_summary_metrics = CollectReadgroupBamQualityMetrics.alignment_summary_metrics
    File read_group_gc_bias_detail_metrics = CollectReadgroupBamQualityMetrics.gc_bias_detail_metrics
    File read_group_gc_bias_pdf = CollectReadgroupBamQualityMetrics.gc_bias_pdf
    File read_group_gc_bias_summary_metrics = CollectReadgroupBamQualityMetrics.gc_bias_summary_metrics

    File selfSM = CheckContamination.selfSM
    File depthSM = CheckContamination.depthSM

    #File validate_cram_report = ValidateCram.report
    File calculate_read_group_checksum_md5 = CalculateReadGroupChecksum.md5_file
    File validate_aggregated_sam_file_report = ValidateAggregatedSamFile.report

    File agg_alignment_summary_metrics = CollectAggregationMetrics.alignment_summary_metrics
    File agg_bait_bias_detail_metrics = CollectAggregationMetrics.bait_bias_detail_metrics
    File agg_bait_bias_summary_metrics = CollectAggregationMetrics.bait_bias_summary_metrics
    File agg_gc_bias_detail_metrics = CollectAggregationMetrics.gc_bias_detail_metrics
    File agg_gc_bias_pdf = CollectAggregationMetrics.gc_bias_pdf
    File agg_gc_bias_summary_metrics = CollectAggregationMetrics.gc_bias_summary_metrics
    File agg_insert_size_histogram_pdf = CollectAggregationMetrics.insert_size_histogram_pdf
    File agg_insert_size_metrics = CollectAggregationMetrics.insert_size_metrics
    File agg_pre_adapter_detail_metrics = CollectAggregationMetrics.pre_adapter_detail_metrics
    File agg_pre_adapter_summary_metrics = CollectAggregationMetrics.pre_adapter_summary_metrics
    File agg_quality_distribution_pdf = CollectAggregationMetrics.quality_distribution_pdf
    File agg_quality_distribution_metrics = CollectAggregationMetrics.quality_distribution_metrics

    File wgs_metrics = CollectWgsMetrics.metrics
    File raw_wgs_metrics = CollectRawWgsMetrics.metrics

    File gvcf_summary_metrics = CollectGvcfCallingMetrics.summary_metrics
    File gvcf_detail_metrics = CollectGvcfCallingMetrics.detail_metrics

    File duplicate_metrics = MarkDuplicates.duplicate_metrics
    File output_bqsr_reports = GatherBqsrReports.output_bqsr_report

    File output_cram = ConvertToCram.output_cram
    File output_cram_index = ConvertToCram.output_cram_index
    File output_cram_md5 = ConvertToCram.output_cram_md5

    File output_vcf = MergeVCFs.output_vcf
    File output_vcf_index = MergeVCFs.output_vcf_index
  }
}
