import "HISAT2.wdl" as HISAT2
import "Picard.wdl" as Picard
import "RSEM.wdl" as RSEM

workflow SmartSeq2SingleCell {
  meta {
    description: "Process SmartSeq2 scRNA-Seq data, include reads alignment, QC metrics collection, and gene expression quantitication"
  }

  # load annotation
  File gtf_file
  File genome_ref_fasta
  File rrna_intervals
  File gene_ref_flat

  # load index
  File hisat2_ref_index
  File hisat2_ref_trans_index
  File rsem_ref_index

  # ref index name
  String hisat2_ref_name
  String hisat2_ref_trans_name

  # samples
  String stranded
  String sample_name
  String output_name
  File fastq1
  File fastq2

  parameter_meta {
    gtf_file: "Gene annotation file in gtf format"
    genome_ref_fasta: "Genome reference in fasta format"
    rrna_intervals: "rRNA interval file required by Picard"
    gene_ref_flat: "Gene refflat file required by Picard"
    hisat2_ref_index: "HISAT2 reference index file in tarball"
    hisat2_ref_trans_index: "HISAT2 transcriptome index file in tarball"
    rsem_ref_index: "RSEM reference index file in tarball"
    hisat2_ref_name: "HISAT2 reference index name"
    hisat2_ref_trans_name: "HISAT2 transcriptome index file name"
    stranded: "Library strand information example values: FR RF NONE"
    sample_name: "Sample name or Cell ID"
    output_name: "Output name, can include path"
    fastq1: "R1 in paired end reads"
    fastq2: "R2 in paired end reads"
  }

  String quality_control_output_basename = output_name + "_qc"

  call HISAT2.HISAT2PairedEnd {
    input:
      hisat2_ref = hisat2_ref_index,
      fastq1 = fastq1,
      fastq2 = fastq2,
      ref_name = hisat2_ref_name,
      sample_name = sample_name,
      output_basename = quality_control_output_basename
  }

  call Picard.CollectMultipleMetrics {
    input:
      aligned_bam = HISAT2PairedEnd.output_bam,
      genome_ref_fasta = genome_ref_fasta,
      output_basename = quality_control_output_basename
  }

  call Picard.CollectRnaMetrics {
    input:
      aligned_bam = HISAT2PairedEnd.output_bam,
      ref_flat = gene_ref_flat,
      rrna_intervals = rrna_intervals,
      output_basename = quality_control_output_basename,
      stranded = stranded,
  }

  call Picard.CollectDuplicationMetrics {
    input:
      aligned_bam = HISAT2PairedEnd.output_bam,
      output_basename = quality_control_output_basename
  }

  String data_output_basename = output_name + "_rsem"

  call HISAT2.HISAT2RSEM as HISAT2Transcriptome {
    input:
      hisat2_ref = hisat2_ref_trans_index,
      fastq1 = fastq1,
      fastq2 = fastq2,
      ref_name = hisat2_ref_trans_name,
      sample_name = sample_name,
      output_basename = data_output_basename,
  }

  call RSEM.RSEMExpression {
    input:
      trans_aligned_bam = HISAT2Transcriptome.output_bam,
      rsem_genome = rsem_ref_index,
      output_basename = data_output_basename,
  }

  output {
    # quality control outputs
    File aligned_bam = HISAT2PairedEnd.output_bam
    File hisat2_met_file = HISAT2PairedEnd.met_file
    File hisat2_log_file = HISAT2PairedEnd.log_file
    File alignment_summary_metrics = CollectMultipleMetrics.alignment_summary_metrics
    File base_call_dist_metrics = CollectMultipleMetrics.base_call_dist_metrics
    File base_call_pdf = CollectMultipleMetrics.base_call_pdf
    File gc_bias_detail_metrics = CollectMultipleMetrics.gc_bias_detail_metrics
    File gc_bias_dist_pdf = CollectMultipleMetrics.gc_bias_dist_pdf
    File gc_bias_summary_metrics = CollectMultipleMetrics.gc_bias_summary_metrics
    File insert_size_hist = CollectMultipleMetrics.insert_size_hist
    File insert_size_metrics = CollectMultipleMetrics.insert_size_metrics
    File quality_distribution_metrics = CollectMultipleMetrics.quality_distribution_metrics
    File quality_distribution_dist_pdf = CollectMultipleMetrics.quality_distribution_dist_pdf
    File quality_by_cycle_metrics = CollectMultipleMetrics.quality_by_cycle_metrics
    File quality_by_cycle_pdf = CollectMultipleMetrics.quality_by_cycle_pdf
    File pre_adapter_details_metrics = CollectMultipleMetrics.pre_adapter_details_metrics
    File bait_bias_detail_metrics = CollectMultipleMetrics.bait_bias_detail_metrics
    File bait_bias_summary_metrics = CollectMultipleMetrics.bait_bias_summary_metrics
    File error_summary_metrics = CollectMultipleMetrics.error_summary_metrics
    File rna_metrics = CollectRnaMetrics.rna_metrics
    File rna_coverage = CollectRnaMetrics.rna_coverage_pdf
    File dedup_metrics = CollectDuplicationMetrics.dedup_metrics

    # data outputs
    File aligned_transcriptome_bam = HISAT2Transcriptome.output_bam
    File hisat2_transcriptome_met_file = HISAT2Transcriptome.met_file
    File hisat2_transcriptome_log_file = HISAT2Transcriptome.log_file
    File rsem_gene_results = RSEMExpression.rsem_gene
    File rsem_isoform_results = RSEMExpression.rsem_isoform
    File rsem_time_log = RSEMExpression.rsem_time
    File rsem_cnt_log = RSEMExpression.rsem_cnt
    File rsem_model_log = RSEMExpression.rsem_model
    File rsem_theta_log = RSEMExpression.rsem_theta
  }
}
