## Copyright Broad Institute, 2018 
## 
## This WDL pipeline implements data pre-processing and initial variant calling (GVCF 
## generation) according to the GATK Best Practices (June 2016) for germline SNP and 
## Indel discovery in human whole-genome sequencing data. 
## 
## Requirements/expectations : 
## - Human whole-genome pair-end sequencing data in unmapped BAM (uBAM) format 
## - One or more read groups, one per uBAM file, all belonging to a single sample (SM) 
## - Input uBAM files must additionally comply with the following requirements: 
## - - filenames all have the same suffix (we use ".unmapped.bam") 
## - - files must pass validation by ValidateSamFile 
## - - reads are provided in query-sorted order 
## - - all reads must have an RG tag 
## - GVCF output names must end in ".g.vcf.gz" 
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
 
import "https://api.firecloud.org/ga4gh/v1/tools/gatk:alignment/versions/4/plain-WDL/descriptor" as Alignment 
import "https://api.firecloud.org/ga4gh/v1/tools/gatk:split-large-readgroup/versions/4/plain-WDL/descriptor" as SplitRG 
import "https://api.firecloud.org/ga4gh/v1/tools/gatk:quality-control/versions/2/plain-WDL/descriptor" as QC 
import "https://api.firecloud.org/ga4gh/v1/tools/gatk:bam-processing/versions/4/plain-WDL/descriptor" as Processing 
import "https://api.firecloud.org/ga4gh/v1/tools/gatk:germline-variant-discovery/versions/3/plain-WDL/descriptor" as Calling 
import "https://api.firecloud.org/ga4gh/v1/tools/gatk:utilities/versions/2/plain-WDL/descriptor" as Utils 
 
# WORKFLOW DEFINITION 
workflow germline_single_sample_workflow { 
 
  File flowcell_unmapped_bams_fofn 
  Array[File] flowcell_unmapped_bams = read_lines(flowcell_unmapped_bams_fofn) 
 
  File contamination_sites_ud 
  File contamination_sites_bed 
  File contamination_sites_mu 
  File wgs_evaluation_interval_list 
  File wgs_coverage_interval_list 
 
  String sample_name 
  String base_file_name 
  String final_vcf_base_name 
  String unmapped_bam_suffix 
 
  File wgs_calling_interval_list 
  Int haplotype_scatter_count 
  Int break_bands_at_multiples_of 
  Int read_length = 250 
 
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
 
  Int preemptible_tries 
  Int agg_preemptible_tries 
 
  Boolean skip_QC 
  Boolean make_gatk4_single_sample_vcf 
  Boolean use_gatk4_haplotype_caller 
 
  Float cutoff_for_large_rg_in_gb = 20.0 
 
  String bwa_commandline="bwa mem -K 100000000 -p -v 3 -t 16 -Y $bash_ref_fasta" 
 
  String recalibrated_bam_basename = base_file_name + ".aligned.duplicates_marked.recalibrated" 
 
  Int compression_level = 2 
 
  # Get the version of BWA to include in the PG record in the header of the BAM produced 
  # by MergeBamAlignment. 
  call Alignment.GetBwaVersion 
 
  # Align flowcell-level unmapped input bams in parallel 
  scatter (unmapped_bam in flowcell_unmapped_bams) { 
 
    Float unmapped_bam_size = size(unmapped_bam, "GB") 
 
    String unmapped_bam_basename = basename(unmapped_bam, unmapped_bam_suffix) 
 
    if (!skip_QC) { 
      # QC the unmapped BAM 
      call QC.CollectQualityYieldMetrics as CollectQualityYieldMetrics { 
        input: 
          input_bam = unmapped_bam, 
          metrics_filename = unmapped_bam_basename + ".unmapped.quality_yield_metrics", 
          preemptible_tries = preemptible_tries 
      } 
    } 
 
    if (unmapped_bam_size > cutoff_for_large_rg_in_gb) { 
      # Split bam into multiple smaller bams, 
      # map reads to reference and recombine into one bam 
      call SplitRG.split_large_readgroup as SplitRG { 
        input: 
          input_bam = unmapped_bam, 
          bwa_commandline = bwa_commandline, 
          bwa_version = GetBwaVersion.version, 
          output_bam_basename = unmapped_bam_basename + ".aligned.unsorted", 
          ref_fasta = ref_fasta, 
          ref_fasta_index = ref_fasta_index, 
          ref_dict = ref_dict, 
          ref_alt = ref_alt, 
          ref_amb = ref_amb, 
          ref_ann = ref_ann, 
          ref_bwt = ref_bwt, 
          ref_pac = ref_pac, 
          ref_sa = ref_sa, 
          compression_level = compression_level, 
          preemptible_tries = preemptible_tries 
      } 
    } 
 
    if (unmapped_bam_size <= cutoff_for_large_rg_in_gb) { 
      # Map reads to reference 
      call Alignment.SamToFastqAndBwaMemAndMba as SamToFastqAndBwaMemAndMba { 
        input: 
          input_bam = unmapped_bam, 
          bwa_commandline = bwa_commandline, 
          output_bam_basename = unmapped_bam_basename + ".aligned.unsorted", 
          ref_fasta = ref_fasta, 
          ref_fasta_index = ref_fasta_index, 
          ref_dict = ref_dict, 
          ref_alt = ref_alt, 
          ref_bwt = ref_bwt, 
          ref_amb = ref_amb, 
          ref_ann = ref_ann, 
          ref_pac = ref_pac, 
          ref_sa = ref_sa, 
          bwa_version = GetBwaVersion.version, 
          compression_level = compression_level, 
          preemptible_tries = preemptible_tries 
      } 
    } 
 
    File output_aligned_bam = select_first([SamToFastqAndBwaMemAndMba.output_bam, SplitRG.aligned_bam]) 
 
    Float mapped_bam_size = size(output_aligned_bam, "GB") 
 
    if (!skip_QC) { 
      # QC the aligned but unsorted readgroup BAM 
      # no reference as the input here is unsorted, providing a reference would cause an error 
      call QC.CollectUnsortedReadgroupBamQualityMetrics as CollectUnsortedReadgroupBamQualityMetrics { 
        input: 
          input_bam = output_aligned_bam, 
          output_bam_prefix = unmapped_bam_basename + ".readgroup", 
          preemptible_tries = preemptible_tries 
      } 
    } 
  } 
 
  # Sum the read group bam sizes to approximate the aggregated bam size 
  call Utils.SumFloats as SumFloats { 
    input: 
      sizes = mapped_bam_size, 
      preemptible_tries = preemptible_tries 
  } 
 
  # Aggregate aligned+merged flowcell BAM files and mark duplicates 
  # We take advantage of the tool's ability to take multiple BAM inputs and write out a single output 
  # to avoid having to spend time just merging BAM files. 
  call Processing.MarkDuplicates as MarkDuplicates { 
    input: 
      input_bams = output_aligned_bam, 
      output_bam_basename = base_file_name + ".aligned.unsorted.duplicates_marked", 
      metrics_filename = base_file_name + ".duplicate_metrics", 
      total_input_size = SumFloats.total_size, 
      compression_level = compression_level, 
      preemptible_tries = agg_preemptible_tries 
  } 
 
  # Sort aggregated+deduped BAM file 
  call Processing.SortSam as SortSampleBam { 
    input: 
      input_bam = MarkDuplicates.output_bam, 
      output_bam_basename = base_file_name + ".aligned.duplicate_marked.sorted", 
      compression_level = compression_level, 
      preemptible_tries = agg_preemptible_tries 
  } 
 
  # Create list of sequences for scatter-gather parallelization 
  call Utils.CreateSequenceGroupingTSV as CreateSequenceGroupingTSV { 
    input: 
      ref_dict = ref_dict, 
      preemptible_tries = preemptible_tries 
  } 
 
  # Estimate level of cross-sample contamination 
  call Processing.CheckContamination as CheckContamination { 
    input: 
      input_bam = SortSampleBam.output_bam, 
      input_bam_index = SortSampleBam.output_bam_index, 
      contamination_sites_ud = contamination_sites_ud, 
      contamination_sites_bed = contamination_sites_bed, 
      contamination_sites_mu = contamination_sites_mu, 
      ref_fasta = ref_fasta, 
      ref_fasta_index = ref_fasta_index, 
      output_prefix = base_file_name + ".preBqsr", 
      preemptible_tries = agg_preemptible_tries, 
      contamination_underestimation_factor = 0.75 
  } 
 
  # We need disk to localize the sharded input and output due to the scatter for BQSR. 
  # If we take the number we are scattering by and reduce by 3 we will have enough disk space 
  # to account for the fact that the data is not split evenly. 
  Int num_of_bqsr_scatters = length(CreateSequenceGroupingTSV.sequence_grouping) 
  Int potential_bqsr_divisor = num_of_bqsr_scatters - 10 
  Int bqsr_divisor = if potential_bqsr_divisor > 1 then potential_bqsr_divisor else 1 
 
  # Perform Base Quality Score Recalibration (BQSR) on the sorted BAM in parallel 
  scatter (subgroup in CreateSequenceGroupingTSV.sequence_grouping) { 
    # Generate the recalibration model by interval 
    call Processing.BaseRecalibrator as BaseRecalibrator { 
      input: 
        input_bam = SortSampleBam.output_bam, 
        recalibration_report_filename = base_file_name + ".recal_data.csv", 
        sequence_group_interval = subgroup, 
        dbSNP_vcf = dbSNP_vcf, 
        dbSNP_vcf_index = dbSNP_vcf_index, 
        known_indels_sites_VCFs = known_indels_sites_VCFs, 
        known_indels_sites_indices = known_indels_sites_indices, 
        ref_dict = ref_dict, 
        ref_fasta = ref_fasta, 
        ref_fasta_index = ref_fasta_index, 
        bqsr_scatter = bqsr_divisor, 
        preemptible_tries = agg_preemptible_tries 
    } 
  } 
 
  # Merge the recalibration reports resulting from by-interval recalibration 
  # The reports are always the same size 
  call Processing.GatherBqsrReports as GatherBqsrReports { 
    input: 
      input_bqsr_reports = BaseRecalibrator.recalibration_report, 
      output_report_filename = base_file_name + ".recal_data.csv", 
      preemptible_tries = preemptible_tries 
  } 
 
  scatter (subgroup in CreateSequenceGroupingTSV.sequence_grouping_with_unmapped) { 
    # Apply the recalibration model by interval 
    call Processing.ApplyBQSR as ApplyBQSR { 
      input: 
        input_bam = SortSampleBam.output_bam, 
        output_bam_basename = recalibrated_bam_basename, 
        recalibration_report = GatherBqsrReports.output_bqsr_report, 
        sequence_group_interval = subgroup, 
        ref_dict = ref_dict, 
        ref_fasta = ref_fasta, 
        ref_fasta_index = ref_fasta_index, 
        bqsr_scatter = bqsr_divisor, 
        compression_level = compression_level, 
        preemptible_tries = agg_preemptible_tries 
    } 
  } 
 
  Float agg_bam_size = size(SortSampleBam.output_bam, "GB") 
 
  # Merge the recalibrated BAM files resulting from by-interval recalibration 
  call Processing.GatherSortedBamFiles as GatherBamFiles { 
    input: 
      input_bams = ApplyBQSR.recalibrated_bam, 
      output_bam_basename = base_file_name, 
      total_input_size = agg_bam_size, 
      compression_level = compression_level, 
      preemptible_tries = agg_preemptible_tries 
  } 
 
  #BQSR bins the qualities which makes a significantly smaller bam 
  Float binned_qual_bam_size = size(GatherBamFiles.output_bam, "GB") 
 
  File? chimerism_metrics 
 
  if (!skip_QC) { 
    # QC the final BAM (consolidated after scattered BQSR) 
    call QC.CollectReadgroupBamQualityMetrics as CollectReadgroupBamQualityMetrics { 
      input: 
        input_bam = GatherBamFiles.output_bam, 
        input_bam_index = GatherBamFiles.output_bam_index, 
        output_bam_prefix = base_file_name + ".readgroup", 
        ref_dict = ref_dict, 
        ref_fasta = ref_fasta, 
        ref_fasta_index = ref_fasta_index, 
        preemptible_tries = agg_preemptible_tries 
    } 
 
    # QC the final BAM some more (no such thing as too much QC) 
    call QC.CollectAggregationMetrics as CollectAggregationMetrics { 
      input: 
        input_bam = GatherBamFiles.output_bam, 
        input_bam_index = GatherBamFiles.output_bam_index, 
        output_bam_prefix = base_file_name, 
        ref_dict = ref_dict, 
        ref_fasta = ref_fasta, 
        ref_fasta_index = ref_fasta_index, 
        preemptible_tries = agg_preemptible_tries 
    } 
 
    # QC the sample WGS metrics (stringent thresholds) 
    call QC.CollectWgsMetrics as CollectWgsMetrics { 
      input: 
        input_bam = GatherBamFiles.output_bam, 
        input_bam_index = GatherBamFiles.output_bam_index, 
        metrics_filename = base_file_name + ".wgs_metrics", 
        ref_fasta = ref_fasta, 
        ref_fasta_index = ref_fasta_index, 
        wgs_coverage_interval_list = wgs_coverage_interval_list, 
        read_length = read_length, 
        preemptible_tries = agg_preemptible_tries 
    } 
 
    # QC the sample raw WGS metrics (common thresholds) 
    call QC.CollectRawWgsMetrics as CollectRawWgsMetrics { 
      input: 
        input_bam = GatherBamFiles.output_bam, 
        input_bam_index = GatherBamFiles.output_bam_index, 
        metrics_filename = base_file_name + ".raw_wgs_metrics", 
        ref_fasta = ref_fasta, 
        ref_fasta_index = ref_fasta_index, 
        wgs_coverage_interval_list = wgs_coverage_interval_list, 
        read_length = read_length, 
        preemptible_tries = agg_preemptible_tries 
    } 
 
    # Generate a checksum per readgroup in the final BAM 
    call QC.CalculateReadGroupChecksum as CalculateReadGroupChecksum { 
      input: 
        input_bam = GatherBamFiles.output_bam, 
        input_bam_index = GatherBamFiles.output_bam_index, 
        read_group_md5_filename = recalibrated_bam_basename + ".bam.read_group_md5", 
        preemptible_tries = agg_preemptible_tries 
    } 
  } 
 
  # ValidateSamFile runs out of memory in mate validation on crazy edge case data, so we want to skip the mate validation 
  # in those cases.  These values set the thresholds for what is considered outside the normal realm of "reasonable" data. 
  Float max_duplication_in_reasonable_sample = 0.30 
  Float max_chimerism_in_reasonable_sample = 0.15 
 
  # Convert the final merged recalibrated BAM file to CRAM format 
  call Utils.ConvertToCram as ConvertToCram { 
    input: 
      input_bam = GatherBamFiles.output_bam, 
      ref_fasta = ref_fasta, 
      ref_fasta_index = ref_fasta_index, 
      output_basename = base_file_name, 
      preemptible_tries = agg_preemptible_tries 
  } 
 
  if (!skip_QC) { 
    # Check whether the data has massively high duplication or chimerism rates 
    call QC.CheckPreValidation as CheckPreValidation { 
      input: 
        duplication_metrics = MarkDuplicates.duplicate_metrics, 
        chimerism_metrics = select_first([CollectAggregationMetrics.alignment_summary_metrics, ""]), 
        max_duplication_in_reasonable_sample = max_duplication_in_reasonable_sample, 
        max_chimerism_in_reasonable_sample = max_chimerism_in_reasonable_sample, 
        preemptible_tries = agg_preemptible_tries 
    } 
  } 
 
  Boolean is_outlier_data = select_first([CheckPreValidation.is_outlier_data, false]) 
 
  # Validate the CRAM file 
  call QC.ValidateSamFile as ValidateCram { 
    input: 
      input_bam = ConvertToCram.output_cram, 
      input_bam_index = ConvertToCram.output_cram_index, 
      report_filename = base_file_name + ".cram.validation_report", 
      ref_dict = ref_dict, 
      ref_fasta = ref_fasta, 
      ref_fasta_index = ref_fasta_index, 
      ignore = ["MISSING_TAG_NM"], 
      max_output = 1000000000, 
      is_outlier_data = is_outlier_data, 
      preemptible_tries = agg_preemptible_tries 
  } 
 
  # Break the calling interval_list into sub-intervals 
  # Perform variant calling on the sub-intervals, and then gather the results 
  call Utils.ScatterIntervalList as ScatterIntervalList { 
    input: 
      interval_list = wgs_calling_interval_list, 
      scatter_count = haplotype_scatter_count, 
      break_bands_at_multiples_of = break_bands_at_multiples_of 
  } 
 
  # We need disk to localize the sharded input and output due to the scatter for HaplotypeCaller. 
  # If we take the number we are scattering by and reduce by 20 we will have enough disk space 
  # to account for the fact that the data is quite uneven across the shards. 
  Int potential_hc_divisor = ScatterIntervalList.interval_count - 20 
  Int hc_divisor = if potential_hc_divisor > 1 then potential_hc_divisor else 1 
 
  # Call variants in parallel over WGS calling intervals 
  scatter (index in range(ScatterIntervalList.interval_count)) { 
 
    if (make_gatk4_single_sample_vcf || use_gatk4_haplotype_caller) { 
      call Calling.HaplotypeCaller_GATK4_VCF as HaplotypeCaller4 { 
        input: 
          contamination = CheckContamination.contamination, 
          input_bam = GatherBamFiles.output_bam, 
          interval_list = ScatterIntervalList.out[index], 
          make_gvcf = !make_gatk4_single_sample_vcf, 
          vcf_basename = base_file_name, 
          ref_dict = ref_dict, 
          ref_fasta = ref_fasta, 
          ref_fasta_index = ref_fasta_index, 
          hc_scatter = hc_divisor, 
          preemptible_tries = agg_preemptible_tries 
      } 
 
      if (make_gatk4_single_sample_vcf) { 
        call Calling.HardFilterVcf as FilterVcf { 
          input: 
            input_vcf = HaplotypeCaller4.output_vcf, 
            input_vcf_index = HaplotypeCaller4.output_vcf_index, 
            vcf_basename = base_file_name, 
            interval_list = ScatterIntervalList.out[index], 
            preemptible_tries = preemptible_tries 
        } 
      } 
    } 
    if (!make_gatk4_single_sample_vcf && !use_gatk4_haplotype_caller) { 
      call Calling.HaplotypeCaller_GATK35_GVCF as HaplotypeCaller3 { 
        input: 
          contamination = CheckContamination.contamination, 
          input_bam = GatherBamFiles.output_bam, 
          interval_list = ScatterIntervalList.out[index], 
          gvcf_basename = base_file_name, 
          ref_dict = ref_dict, 
          ref_fasta = ref_fasta, 
          ref_fasta_index = ref_fasta_index, 
          hc_scatter = hc_divisor, 
          preemptible_tries = agg_preemptible_tries 
       } 
     } 
 
     File merge_input = select_first([FilterVcf.output_vcf, HaplotypeCaller4.output_vcf, HaplotypeCaller3.output_gvcf]) 
     File merge_input_index = select_first([FilterVcf.output_vcf_index, HaplotypeCaller4.output_vcf_index, HaplotypeCaller3.output_gvcf_index]) 
  } 
 
  String name_token = if make_gatk4_single_sample_vcf then ".filtered" else ".g" 
 
  # Combine by-interval VCFs into a single sample VCF file 
  call Calling.MergeVCFs as MergeVCFs { 
    input: 
      input_vcfs = merge_input, 
      input_vcfs_indexes = merge_input_index, 
      output_vcf_name = final_vcf_base_name + name_token + ".vcf.gz", 
      preemptible_tries = agg_preemptible_tries 
  } 
 
  # Outputs that will be retained when execution is complete 
  output { 
    Array[File?] quality_yield_metrics = CollectQualityYieldMetrics.quality_yield_metrics 
 
    Array[File?] unsorted_read_group_base_distribution_by_cycle_pdf = CollectUnsortedReadgroupBamQualityMetrics.base_distribution_by_cycle_pdf 
    Array[File?] unsorted_read_group_base_distribution_by_cycle_metrics = CollectUnsortedReadgroupBamQualityMetrics.base_distribution_by_cycle_metrics 
    Array[File?] unsorted_read_group_insert_size_histogram_pdf = CollectUnsortedReadgroupBamQualityMetrics.insert_size_histogram_pdf 
    Array[File?] unsorted_read_group_insert_size_metrics = CollectUnsortedReadgroupBamQualityMetrics.insert_size_metrics 
    Array[File?] unsorted_read_group_quality_by_cycle_pdf = CollectUnsortedReadgroupBamQualityMetrics.quality_by_cycle_pdf 
    Array[File?] unsorted_read_group_quality_by_cycle_metrics = CollectUnsortedReadgroupBamQualityMetrics.quality_by_cycle_metrics 
    Array[File?] unsorted_read_group_quality_distribution_pdf = CollectUnsortedReadgroupBamQualityMetrics.quality_distribution_pdf 
    Array[File?] unsorted_read_group_quality_distribution_metrics = CollectUnsortedReadgroupBamQualityMetrics.quality_distribution_metrics 
 
    File? read_group_alignment_summary_metrics = CollectReadgroupBamQualityMetrics.alignment_summary_metrics 
    File? read_group_gc_bias_detail_metrics = CollectReadgroupBamQualityMetrics.gc_bias_detail_metrics 
    File? read_group_gc_bias_pdf = CollectReadgroupBamQualityMetrics.gc_bias_pdf 
    File? read_group_gc_bias_summary_metrics = CollectReadgroupBamQualityMetrics.gc_bias_summary_metrics 
 
    File selfSM = CheckContamination.selfSM 
    Float contamination = CheckContamination.contamination 
 
    File? calculate_read_group_checksum_md5 = CalculateReadGroupChecksum.md5_file 
 
    File? agg_alignment_summary_metrics = CollectAggregationMetrics.alignment_summary_metrics 
    File? agg_bait_bias_detail_metrics = CollectAggregationMetrics.bait_bias_detail_metrics 
    File? agg_bait_bias_summary_metrics = CollectAggregationMetrics.bait_bias_summary_metrics 
    File? agg_gc_bias_detail_metrics = CollectAggregationMetrics.gc_bias_detail_metrics 
    File? agg_gc_bias_pdf = CollectAggregationMetrics.gc_bias_pdf 
    File? agg_gc_bias_summary_metrics = CollectAggregationMetrics.gc_bias_summary_metrics 
    File? agg_insert_size_histogram_pdf = CollectAggregationMetrics.insert_size_histogram_pdf 
    File? agg_insert_size_metrics = CollectAggregationMetrics.insert_size_metrics 
    File? agg_pre_adapter_detail_metrics = CollectAggregationMetrics.pre_adapter_detail_metrics 
    File? agg_pre_adapter_summary_metrics = CollectAggregationMetrics.pre_adapter_summary_metrics 
    File? agg_quality_distribution_pdf = CollectAggregationMetrics.quality_distribution_pdf 
    File? agg_quality_distribution_metrics = CollectAggregationMetrics.quality_distribution_metrics 
 
    File? wgs_metrics = CollectWgsMetrics.metrics 
    File? raw_wgs_metrics = CollectRawWgsMetrics.metrics 
 
    File duplicate_metrics = MarkDuplicates.duplicate_metrics 
    File output_bqsr_reports = GatherBqsrReports.output_bqsr_report 
 
    File output_cram = ConvertToCram.output_cram 
    File output_cram_index = ConvertToCram.output_cram_index 
    File output_cram_md5 = ConvertToCram.output_cram_md5 
 
    File validate_cram_file_report = ValidateCram.report 
 
    File output_vcf = MergeVCFs.output_vcf 
    File output_vcf_index = MergeVCFs.output_vcf_index 
  } 
}
