task StripSuffix {
  String filename
  String suffix

  # Use python to strip the leading path and trailing extension
  # from the input file name
  # /path/to/foo.txt -> foo
  command <<<
    python <<CODE
    from os.path import basename, splitext
    out = basename("${filename}").rstrip("${suffix}")
    print out
    CODE
  >>>
  runtime {
    docker: "python:2.7"
    memory: "2 GB"
  }
  output {
    String base_name = read_string(stdout())
  }
}

task CollectQualityYieldMetrics {
  File input_bam
  String metrics_filename

  command {
    java -Xmx128m -jar /usr/gitc/picard.jar \
      CollectQualityYieldMetrics \
      INPUT=${input_bam} \
      OQ=true \
      OUTPUT=${metrics_filename}
  }
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:1.953"
    disks: "local-disk 100 SSD"
    memory: "2 GB"
  }
  output {
    File metrics = "${metrics_filename}"
  }
}

task SamToFastqAndBwaMem {
  File input_bam
  String output_bam_basename
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  File ref_amb
  File ref_ann
  File ref_bwt
  File ref_pac
  File ref_sa

  command {
    set -o pipefail
    java -Xmx3000m -jar /usr/gitc/picard.jar \
      SamToFastq \
      INPUT=${input_bam} \
      FASTQ=/dev/stdout \
      INTERLEAVE=true \
      CLIPPING_ATTRIBUTE=XT \
      CLIPPING_ACTION=2 \
      NON_PF=true |
    /usr/gitc/bwa mem -M -p -v 3 -t 16 ${ref_fasta} /dev/stdin > ${output_bam_basename}.bam
  }
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:1.953"
    memory: "14 GB"
    cpu: "16"
    disks: "local-disk 200 SSD"
  }
  output {
    File output_bam = "${output_bam_basename}.bam"
  }
}

task MergeBamAlignment {
  File unmapped_bam
  File aligned_bam
  String output_bam_basename
  File ref_fasta
  File ref_fasta_index
  File ref_dict  

  command {
    java -Xmx3000m -jar /usr/gitc/picard.jar \
      MergeBamAlignment \
      VALIDATION_STRINGENCY=SILENT \
      CREATE_INDEX=true \
      EXPECTED_ORIENTATIONS=FR \
      ATTRIBUTES_TO_RETAIN=X0 \
      ALIGNED_BAM=${aligned_bam} \
      UNMAPPED_BAM=${unmapped_bam} \
      OUTPUT=${output_bam_basename}.bam \
      REFERENCE_SEQUENCE=${ref_fasta} \
      PAIRED_RUN=true \
      IS_BISULFITE_SEQUENCE=false \
      ALIGNED_READS_ONLY=false \
      CLIP_ADAPTERS=false \
      MAX_RECORDS_IN_RAM=2000000 \
      ADD_MATE_CIGAR=true \
      MAX_INSERTIONS_OR_DELETIONS=-1 \
      PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
      PROGRAM_RECORD_ID="bwamem" \
      PROGRAM_GROUP_VERSION="0.7.7-r441" \
      PROGRAM_GROUP_COMMAND_LINE="bwa mem -M -p -v 3 -t 16 ${ref_fasta}" \
      PROGRAM_GROUP_NAME="bwamem"
  }
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:1.953"
    memory: "3500 MB"
    cpu: "1"
    disks: "local-disk 200 SSD"
  }
  output {
    File output_bam = "${output_bam_basename}.bam"
    File output_bam_index = "${output_bam_basename}.bai"
  }  
}

task CollectMultipleMetrics {
  File input_bam
  File input_bam_index
  String output_bam_prefix
  Array[String] programs
  Array[String] metric_accumulation_levels
  File ref_dict
  File ref_fasta
  File ref_fasta_index

  command {
    java -Xmx5000m -jar /usr/gitc/picard.jar \
      CollectMultipleMetrics \
      INPUT=${input_bam} \
      REFERENCE_SEQUENCE=${ref_fasta} \
      OUTPUT=${output_bam_prefix} \
      ASSUME_SORTED=true \
      PROGRAM=${sep=' PROGRAM=' programs} \
      METRIC_ACCUMULATION_LEVEL=${sep=' METRIC_ACCUMULATION_LEVEL=' metric_accumulation_levels}
  }
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:1.953"
    memory: "7 GB"
    disks: "local-disk 200 SSD"
  }
  output {
    Array[File] metrics = glob("${output_bam_prefix}.*")
  }
}

task CrossCheckFingerprints {
  Array[File] input_bams
  Array[File] input_bam_indexes
  File haplotype_database_file
  String metrics_filename

  command {  
    java -Dsamjdk.buffer_size=131072 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx2000m \
     -jar /usr/gitc/picard.jar \
     CrosscheckReadGroupFingerprints \
     OUTPUT=${metrics_filename} \
     HAPLOTYPE_MAP=${haplotype_database_file} \
     EXPECT_ALL_READ_GROUPS_TO_MATCH=true \
     INPUT=${sep=' INPUT=' input_bams} \
     LOD_THRESHOLD=-20.0
  }
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:1.953"
    memory: "2 GB"
    disks: "local-disk 200 SSD"
  }
  output {
    File metrics = "${metrics_filename}"
  }
}

task CheckReadGroupFingerprints {
  File input_bam
  File input_bam_index
  File haplotype_database_file
  File genotypes
  String metrics_basename
  String sample
  command <<<
  if [ -s ${genotypes} ]; then
    java -Dsamjdk.buffer_size=131072 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx1024m  \
    -jar /usr/gitc/picard.jar \
    CheckFingerprint \
    INPUT=${input_bam} \
    OUTPUT=${metrics_basename} \
    GENOTYPES=${genotypes} \
    HAPLOTYPE_MAP=${haplotype_database_file} \
    SAMPLE_ALIAS=${sample}
  else
    echo "No fingerprint found. Skipping Fingerprint check"
  fi
  >>>
 runtime {
    docker: "broadinstitute/genomes-in-the-cloud:1.953"
    memory: "1 GB"
    disks: "local-disk 200 SSD"
  }
  output {
    File summary_metrics = "${metrics_basename}.fingerprinting_summary_metrics"
    File detail_metrics = "${metrics_basename}.fingerprinting_detail_metrics"
  }
}

task MarkDuplicates {
  Array[File] input_bams
  Array[File] input_bam_indexes
  String output_bam_basename
  String metrics_filename

  command {
    java -Xmx4000m -jar /usr/gitc/picard.jar \
      MarkDuplicates \
      CREATE_INDEX=true \
      CREATE_MD5_FILE=true \
      INPUT=${sep=' INPUT=' input_bams} \
      OUTPUT=${output_bam_basename}.bam \
      METRICS_FILE=${metrics_filename} \
      VALIDATION_STRINGENCY=SILENT \
      OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500
  }
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:1.953"
    memory: "7 GB"
    disks: "local-disk 300 SSD"
  }
  output {
    File output_bam = "${output_bam_basename}.bam"
    File output_bam_index = "${output_bam_basename}.bai"
    File output_bam_checksum = "${output_bam_basename}.bam.md5"
    File duplicate_metrics = "${metrics_filename}"
  }
}

task CreateSequenceGroupingTSV {
  File ref_dict

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
                sequence_tuple_list.append((line_split[1].split(":")[1], int(line_split[2].split(":")[1])))
        longest_sequence = sorted(sequence_tuple_list, key=lambda x: x[1], reverse=True)[0][1]

    # initialize the tsv string with the first sequence
    tsv_string = sequence_tuple_list[0][0]
    temp_size = sequence_tuple_list[0][1]
    for sequence_tuple in sequence_tuple_list[1:]:
        if temp_size + sequence_tuple[1] <= longest_sequence:
            temp_size += sequence_tuple[1]
            tsv_string += "\t" + sequence_tuple[0]
        else:
            tsv_string += "\n" + sequence_tuple[0]
            temp_size = sequence_tuple[1]

    print tsv_string
    CODE
  >>>
  runtime {
    docker: "python:2.7"
    memory: "2 GB"
  }
  output {
    Array[Array[String]] sequence_grouping = read_tsv(stdout())
  }
}

task BaseRecalibrator {
  File input_bam
  File input_bam_index
  String recalibration_report_filename
  Array[String] sequence_group_interval
  File dbSNP_vcf
  File dbSNP_vcf_index
  File known_snps_sites_vcf
  File known_snps_sites_vcf_index
  File known_indels_sites_vcf
  File known_indels_sites_vcf_index
  File ref_dict
  File ref_fasta
  File ref_fasta_index

  command {
    java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx4000m \
      -jar /usr/gitc/GenomeAnalysisTK-3.4-g3c929b0.jar \
      -T BaseRecalibrator \
      --disable_auto_index_creation_and_locking_when_reading_rods \
      -U \
      -R ${ref_fasta} \
      -I ${input_bam} \
      --useOriginalQualities \
      -o ${recalibration_report_filename} \
      -knownSites ${dbSNP_vcf} \
      -knownSites ${known_snps_sites_vcf} \
      -knownSites ${known_indels_sites_vcf} \
      -L ${sep=" -L " sequence_group_interval}
  }
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:1.953"
    memory: "6 GB"
    disks: "local-disk 200 SSD"
  }
  output {
    File recalibration_report = "${recalibration_report_filename}"
  }
}

task BqsrPrintReads {
  File input_bam
  File input_bam_index
  String output_bam_basename
  File recalibration_report
  Array[String] sequence_group_interval
  File ref_dict
  File ref_fasta
  File ref_fasta_index

  command {
    java -Dsamjdk.use_async_io=true -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx3000m \
      -jar /usr/gitc/GenomeAnalysisTK-3.4-g3c929b0.jar \
      -T PrintReads \
      --disable_auto_index_creation_and_locking_when_reading_rods \
      --generate_md5 \
      -U \
      -R ${ref_fasta} \
      -I ${input_bam} \
      --useOriginalQualities \
      -o ${output_bam_basename}.bam \
      --disable_indel_quals \
      -BQSR ${recalibration_report} \
      --emit_original_quals \
      -L ${sep=" -L " sequence_group_interval}
  }
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:1.953"
    memory: "3500 MB"
    disks: "local-disk 200 SSD"
  }
  output {
    File recalibrated_bam = "${output_bam_basename}.bam"
    File recalibrated_bam_checksum = "${output_bam_basename}.bam.md5"
  }
}

task GatherBqsrReports {
  Array[File] input_bqsr_reports
  String output_report_filename

  command {
    java -Xmx3000m -cp /usr/gitc/GatherBqsrReports-GATK.jar \
      org.broadinstitute.gatk.tools.GatherBqsrReports \
      INPUT=${sep=' INPUT=' input_bqsr_reports} \
      OUTPUT=${output_report_filename}
    }
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:GatherBqsrReport"
    memory: "3500 MB"
    disks: "local-disk 100 SSD"
  }
  output {
    File output_bqsr_report = "${output_report_filename}"
  }
}

task GatherBamFiles {
  Array[File] input_bams
  File input_unmapped_reads_bam
  String output_bam_basename

  command {
    java -Xmx4000m -jar /usr/gitc/picard.jar \
      GatherBamFiles \
      INPUT=${sep=' INPUT=' input_bams} \
      INPUT=${input_unmapped_reads_bam} \
      OUTPUT=${output_bam_basename}.bam \
      CREATE_INDEX=true \
      CREATE_MD5_FILE=true

    }
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:1.953"
    memory: "7 GB"
    disks: "local-disk 400 SSD"
  }
  output {
    File output_bam = "${output_bam_basename}.bam"
    File output_bam_index = "${output_bam_basename}.bai"
    File output_bam_md5 = "${output_bam_basename}.bam.md5"
  }
}

task ValidateSamFile {
  File input_bam
  File input_bam_index
  String report_filename
  File ref_dict
  File ref_fasta
  File ref_fasta_index

  command {
    java -Xmx4000m -jar /usr/gitc/picard.jar \
      ValidateSamFile \
      INPUT=${input_bam} \
      OUTPUT=${report_filename} \
      REFERENCE_SEQUENCE=${ref_fasta} \
      MODE=SUMMARY \
      IS_BISULFITE_SEQUENCED=false \
      CREATE_MD5_FILE=false
  }
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:1.953"
    memory: "7 GB"
    disks: "local-disk 200 SSD"
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
  File ref_fasta
  File ref_fasta_index

  command {
    java -Xmx4000m -jar /usr/gitc/picard.jar \
      CollectWgsMetrics \
      INPUT=${input_bam} \
      VALIDATION_STRINGENCY=SILENT \
      REFERENCE_SEQUENCE=${ref_fasta} \
      OUTPUT=${metrics_filename}
  }
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:1.953"
    memory: "7 GB"
    disks: "local-disk 200 SSD"
  }
  output {
    File metrics = "${metrics_filename}"
  }
}

task CollectRawWgsMetrics {
  File input_bam
  File input_bam_index
  String metrics_filename
  File ref_fasta
  File ref_fasta_index

  command {
    java -Xmx4000m -jar /usr/gitc/picard.jar \
      CollectRawWgsMetrics \
      INPUT=${input_bam} \
      VALIDATION_STRINGENCY=SILENT \
      REFERENCE_SEQUENCE=${ref_fasta} \
      OUTPUT=${metrics_filename}
  }
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:1.953"
    memory: "7 GB"
    disks: "local-disk 200 SSD"
  }
  output {
    File metrics = "${metrics_filename}"
  }
}

task CalculateReadGroupChecksum {
  File input_bam
  File input_bam_index
  String read_group_md5_filename

  command {
    java -Xmx4000m -jar /usr/gitc/picard.jar \
      CalculateReadGroupChecksum \
      INPUT=${input_bam} \
      OUTPUT=${read_group_md5_filename}
  }
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:1.953"
    memory: "7 GB"
    disks: "local-disk 200 SSD"
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

  # Having to do this as a 2-step command in heredoc syntax, adding a python step to read the metrics
  # This is a hack until read_object() is supported by Cromwell.
  # It relies on knowing that there is only one data row in the 2-row selfSM TSV file
  # Piping output of verifyBamId to /dev/null so only stdout is from the python command
  command <<<
    /usr/gitc/verifyBamID \
    --verbose \
    --ignoreRG \
    --vcf ${contamination_sites_vcf} \
    --out ${output_prefix} \
    --bam ${input_bam} \
    1>/dev/null

    python3 <<CODE
    import csv
    with open('${output_prefix}.selfSM') as selfSM:
      reader = csv.DictReader(selfSM, delimiter='\t')
      for row in reader:
        print(float(row["FREEMIX"])/0.75)
    CODE
  >>>

  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:1.953"
    memory: "2 GB"
    disks: "local-disk 200 SSD"
  }
  output {
    File selfSM = "${output_prefix}.selfSM"
    File depthSM = "${output_prefix}.depthSM"
    File log = "${output_prefix}.log"

    # I would like to do the following, however:
    # The objct is read as a string
    # explicit string->float coercion via float(), as shown below, is supported by Cromwell
    # the interim value cannot be store as a string and then assigned to a float. Variables intialized in output cannot be dereferenced in output.
    # Float contamination = float(read_object(${output_prefix} + ".selfSM").FREEMIX) / 0.75

    # In the interim, get the value from the python hack above:
    Float contamination = read_float(stdout())
  }
}

task ScatterIntervalList {
  File interval_list
  Int scatter_count
  Int break_bands_at_multiples_of

  command <<<
    mkdir out
    java -Xmx1g -jar /usr/gitc/picard.jar \
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
    for i, interval in enumerate(glob.glob("out/*/*.interval_list")):
      (directory, filename) = os.path.split(interval)
      newName = os.path.join(directory, str(i) + filename)
      os.rename(interval, newName)
    CODE
  >>>
  output {
    Array[File] out = glob("out/*/*.interval_list")
  }

  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:1.953"
    memory: "2 GB"
  }
}

task HaplotypeCaller {
  File input_bam
  File input_bam_index
  Float contamination
  File interval_list
  String gvcf_basename
  File ref_dict
  File ref_fasta
  File ref_fasta_index

  # tried to find lowest memory variable where it would still work, might change once tested on JES
  command {
    java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx6000m \
      -jar /usr/gitc/GenomeAnalysisTK-3.4-g3c929b0.jar \
      -T HaplotypeCaller \
      -R ${ref_fasta} \
      -o ${gvcf_basename}.vcf.gz \
      -I ${input_bam} \
      -L ${interval_list} \
      -ERC GVCF \
      --max_alternate_alleles 3 \
      -variant_index_parameter 128000 \
      -variant_index_type LINEAR \
      -contamination ${contamination} \
      --read_filter OverclippedRead
  }
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:1.953"
    memory: "8 GB"
    cpu: "1"
    disks: "local-disk 200 SSD"
  }
  output {
    File output_gvcf = "${gvcf_basename}.vcf.gz"
    File output_gvcf_index = "${gvcf_basename}.vcf.gz.tbi"
  }
}

task GatherVCFs {
  Array[File] input_vcfs
  Array[File] input_vcfs_indexes
  String output_vcf_basename

  # using MergeVcfs instead of GatherVcfs so we can create indicies
  # WARNING	2015-10-28 15:01:48	GatherVcfs	Index creation not currently supported when gathering block compressed VCFs.
  command {
    java -Xmx2g -jar /usr/gitc/picard.jar \
    MergeVcfs \
    INPUT=${sep=' INPUT=' input_vcfs} \
    OUTPUT="${output_vcf_basename}.vcf.gz"
  }
  output {
    File output_vcf = "${output_vcf_basename}.vcf.gz"
    File output_vcf_index = "${output_vcf_basename}.vcf.gz.tbi"
  }
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:1.953"
    memory: "3 GB"
    disks: "local-disk 200 SSD"
  }
}

workflow PairedEndSingleSampleWorkflow {
  String sample_name
  Array[File] flowcell_unmapped_bams
  String unmapped_bam_suffix
  File calling_interval_list
  File dbSNP_vcf
  File dbSNP_vcf_index
  File known_snps_sites_vcf
  File known_snps_sites_vcf_index
  File known_indels_sites_vcf
  File known_indels_sites_vcf_index
  File contamination_sites_vcf
  File contamination_sites_vcf_index
  File fingerprint_genotypes_file # if this file is empty (0-length) the workflow should not do fingerprint comparison (as there are no fingerprints for the sample)
  File haplotype_database_file
  Int haplotype_scatter_count
  Int break_bands_at_multiples_of
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  File ref_bwt
  File ref_sa
  File ref_amb
  File ref_ann
  File ref_pac

  String recalibrated_bam_basename = sample_name + ".aligned.duplicates_marked.recalibrated"

  # Align flowcell-level unmapped input bams in parallel
  scatter (unmapped_bam in flowcell_unmapped_bams) {
    call StripSuffix as StripBamExtension {
      input:
        filename = unmapped_bam,
        suffix = ".bam"
    }
    call CollectQualityYieldMetrics {
      input:
        input_bam = unmapped_bam,
        metrics_filename = StripBamExtension.base_name + ".quality_yield_metrics"

    }
    call StripSuffix as StripUnmappedBamSuffix {
      input:
        filename = unmapped_bam,
        suffix = unmapped_bam_suffix
    }
    call SamToFastqAndBwaMem {
      input:
        input_bam = unmapped_bam,
        output_bam_basename = StripUnmappedBamSuffix.base_name + ".unmerged",
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        ref_bwt = ref_bwt,
        ref_amb = ref_amb,
        ref_ann = ref_ann,
        ref_pac = ref_pac,
        ref_sa = ref_sa
    }
    call MergeBamAlignment {
      input:      
        unmapped_bam = unmapped_bam,
        aligned_bam = SamToFastqAndBwaMem.output_bam,
        output_bam_basename = StripUnmappedBamSuffix.base_name + ".aligned",
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict        
    }   

    # called to help in finding problems early. 
    # if too time consuming and not helpful, can be removed. 
    call ValidateSamFile as ValidateReadGroupSamFile {
    input:
      input_bam = ReadGroupMarkDuplicates.output_bam,
      input_bam_index = ReadGroupMarkDuplicates.output_bam_index,
      report_filename = StripUnmappedBamSuffix.base_name + ".validation_report",
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index
    }
 
    # Need to markduplicates here so that the read-group metrics will not be affected by duplicates. 
    # Also it is good to have read-group duplication metrics for QC. Once we are convinced that the RapidQC 
    # metrics agree with these metrics even in the worst of times, the following two step can be removed.
    call MarkDuplicates as ReadGroupMarkDuplicates {
      input:
        input_bams = MergeBamAlignment.output_bam,
        input_bam_indexes = MergeBamAlignment.output_bam_index,
        output_bam_basename = StripUnmappedBamSuffix.base_name + ".aligned.duplicates_marked",
        metrics_filename = StripUnmappedBamSuffix.base_name + ".duplicate_metrics"
    }

    call CheckReadGroupFingerprints {
    input:
      input_bam = ReadGroupMarkDuplicates.output_bam,
      input_bam_index = ReadGroupMarkDuplicates.output_bam_index,
      haplotype_database_file = haplotype_database_file,
      genotypes = fingerprint_genotypes_file,
      metrics_basename = StripUnmappedBamSuffix.base_name,
      sample = sample_name
    }

    call CollectMultipleMetrics as CollectReadgroupBamQualityMetrics {
    input:
      input_bam = ReadGroupMarkDuplicates.output_bam,
      input_bam_index = ReadGroupMarkDuplicates.output_bam_index,
      output_bam_prefix = StripUnmappedBamSuffix.base_name,
      programs = ["null",                                    # 'null' clears the default list
                  "CollectAlignmentSummaryMetrics",
                  "CollectBaseDistributionByCycle",
                  "CollectInsertSizeMetrics",
                  "MeanQualityByCycle",
                  "QualityScoreDistribution",
                  "CollectGcBiasMetrics"],
      metric_accumulation_levels = ["null", "ALL_READS"],
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index
    }
  }

  call CrossCheckFingerprints {
  input:
    input_bams = ReadGroupMarkDuplicates.output_bam,
    input_bam_indexes = ReadGroupMarkDuplicates.output_bam_index,
    haplotype_database_file = haplotype_database_file,
    metrics_filename = sample_name + ".crosscheck"
  }

  # Aggregate aligned+merged flowcell bams, mark duplicates, and recalibrate qualities
  # This can use the non-duplicate-Marked readgroup bams since it will mark duplicates again.
  call MarkDuplicates as AggregatedBamMarkDuplicates {
    input:
      input_bams = MergeBamAlignment.output_bam,
      input_bam_indexes = MergeBamAlignment.output_bam_index,
      output_bam_basename = sample_name + ".aligned.duplicates_marked",
      metrics_filename = sample_name + ".duplicate_metrics"
  }
  call CreateSequenceGroupingTSV {
    input:
      ref_dict = ref_dict
  }
  call CheckContamination as PreBqsrCheckContamination {
    input:
      input_bam = AggregatedBamMarkDuplicates.output_bam,
      input_bam_index = AggregatedBamMarkDuplicates.output_bam_index,
      contamination_sites_vcf = contamination_sites_vcf,
      contamination_sites_vcf_index = contamination_sites_vcf_index,
      output_prefix = sample_name + ".preBqsr"
  }

  scatter (subgroup in CreateSequenceGroupingTSV.sequence_grouping) {
    call BaseRecalibrator {
      input:
        input_bam = AggregatedBamMarkDuplicates.output_bam,
        input_bam_index = AggregatedBamMarkDuplicates.output_bam_index,
        recalibration_report_filename = sample_name + ".recal_data.csv",
        sequence_group_interval = subgroup,
        dbSNP_vcf = dbSNP_vcf,
        dbSNP_vcf_index = dbSNP_vcf_index,
        known_snps_sites_vcf = known_snps_sites_vcf,
        known_snps_sites_vcf_index = known_snps_sites_vcf_index,
        known_indels_sites_vcf = known_indels_sites_vcf,
        known_indels_sites_vcf_index = known_indels_sites_vcf_index,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index
    }
    call BqsrPrintReads as ScatteredPrintReads {
      input:
        input_bam = AggregatedBamMarkDuplicates.output_bam,
        input_bam_index = AggregatedBamMarkDuplicates.output_bam_index,
        output_bam_basename = recalibrated_bam_basename,
        recalibration_report = GatherBqsrReports.output_bqsr_report,
        sequence_group_interval = subgroup,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index
    }
  }
  call GatherBqsrReports{
    input:
      input_bqsr_reports = BaseRecalibrator.recalibration_report,
      output_report_filename = sample_name + ".recal_data.csv"
  }
  # this printreads uses "-L unmapped" as its interval to grab all paired unmapped reads so they can be added to the final bam
  Array[String] unmapped_group_interval = ["unmapped"]
  call BqsrPrintReads as PrintUnmappedReads {
    input:
      input_bam = AggregatedBamMarkDuplicates.output_bam,
      input_bam_index = AggregatedBamMarkDuplicates.output_bam_index,
      output_bam_basename = recalibrated_bam_basename,
      recalibration_report = GatherBqsrReports.output_bqsr_report,
      sequence_group_interval = unmapped_group_interval,
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index
  }

  ## TODO when capability of adding elements to arrays, can just have one array as an input and add the output of the above task to the scattered printreads bams
  call GatherBamFiles {
    input:
      input_bams = ScatteredPrintReads.recalibrated_bam,
      input_unmapped_reads_bam = PrintUnmappedReads.recalibrated_bam,
      output_bam_basename = sample_name
  }
  call ValidateSamFile as ValidateAggregatedSamFile{
    input:
      input_bam = GatherBamFiles.output_bam,
      input_bam_index = GatherBamFiles.output_bam_index,
      report_filename = sample_name + ".validation_report",
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index
  }
  call CollectMultipleMetrics as CollectAggregationMetrics {
    input:
      input_bam = GatherBamFiles.output_bam,
      input_bam_index = GatherBamFiles.output_bam_index,
      output_bam_prefix = sample_name,
      programs = ["null",                   # 'null' clears the default list
                  "CollectAlignmentSummaryMetrics",
                  "CollectInsertSizeMetrics",
                  "CollectSequencingArtifactMetrics",
                  "CollectGcBiasMetrics"],
      metric_accumulation_levels = ["null", "SAMPLE", "LIBRARY"],
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index
  }
  call CollectWgsMetrics {
    input:
      input_bam = GatherBamFiles.output_bam,
      input_bam_index = GatherBamFiles.output_bam_index,
      metrics_filename = sample_name + ".wgs_metrics",
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index
  }
  call CollectRawWgsMetrics {
    input:
      input_bam = GatherBamFiles.output_bam,
      input_bam_index = GatherBamFiles.output_bam_index,
      metrics_filename = sample_name + ".raw_wgs_metrics",
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index
  }
  call CalculateReadGroupChecksum {
    input:
      input_bam = GatherBamFiles.output_bam,
      input_bam_index = GatherBamFiles.output_bam_index,
      read_group_md5_filename = recalibrated_bam_basename + ".bam.read_group_md5"
  }
  
  #TODO: once we tie-out the contamination results from pre and post bqsr, this task can be 
  #removed.
  call CheckContamination {
    input:
      input_bam = GatherBamFiles.output_bam,
      input_bam_index = GatherBamFiles.output_bam_index,
      contamination_sites_vcf = contamination_sites_vcf,
      contamination_sites_vcf_index = contamination_sites_vcf_index,
      output_prefix = sample_name
  }

  # Break the calling interval_list into sub-intervals
  # Perform variant calling on the sub-intervals, and then gather the results
  call ScatterIntervalList {
    input:
      interval_list = calling_interval_list,
      scatter_count = haplotype_scatter_count,
      break_bands_at_multiples_of = break_bands_at_multiples_of
  }
  scatter (subInterval in ScatterIntervalList.out) {
     call HaplotypeCaller {
       input:
         input_bam = GatherBamFiles.output_bam,
         input_bam_index = GatherBamFiles.output_bam_index,

         # once we see that Pre- and Post-BQSR check contamination is similar, we can uncomment the next line
         # and remove the one the follows
         # contamination = PreBqsrCheckContamination.contamination,
         contamination = CheckContamination.contamination,
         interval_list = subInterval,
         gvcf_basename = sample_name,
         ref_dict = ref_dict,
         ref_fasta = ref_fasta,
         ref_fasta_index = ref_fasta_index
     }
  }
  call GatherVCFs {
    input:
      input_vcfs = HaplotypeCaller.output_gvcf,
      input_vcfs_indexes = HaplotypeCaller.output_gvcf_index,
      output_vcf_basename = sample_name + ".final"
  }
  output {
    CollectQualityYieldMetrics.*
    ValidateReadGroupSamFile.*
    CollectReadgroupBamQualityMetrics.*
    ReadGroupMarkDuplicates.duplicate_metrics
    CheckReadGroupFingerprints.*
    CrossCheckFingerprints.*
    AggregatedBamMarkDuplicates.duplicate_metrics
    GatherBqsrReports.*
    GatherBamFiles.*
    CalculateReadGroupChecksum.*
    ValidateAggregatedSamFile.*
    CollectAggregationMetrics.*
    CollectWgsMetrics.*
    CollectRawWgsMetrics.*
    CheckContamination.*
    GatherVCFs.*
  }
}
