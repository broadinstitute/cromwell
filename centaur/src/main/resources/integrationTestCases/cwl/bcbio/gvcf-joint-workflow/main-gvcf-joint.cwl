cwlVersion: v1.0
class: Workflow
hints: []
inputs:
- id: config__algorithm__align_split_size
  type:
    items:
    - string
    - 'null'
    - boolean
    type: array
- id: files
  type:
    items:
      items: File
      type: array
    type: array
- id: config__algorithm__trim_reads
  type:
    items:
    - string
    - 'null'
    - boolean
    type: array
- id: reference__fasta__base
  secondaryFiles:
  - .fai
  - ^.dict
  type:
    items: File
    type: array
- id: metadata__phenotype
  type:
    items: string
    type: array
- id: config__algorithm__vcfanno
  type:
    items:
      items: File
      type: array
    type: array
- id: resources
  type:
    items: string
    type: array
- id: config__algorithm__variantcaller
  type:
    items:
      items: string
      type: array
    type: array
- id: config__algorithm__adapters
  type:
    items:
    - 'null'
    - items:
      - 'null'
      - string
      type: array
    type: array
- id: genome_resources__variation__1000g
  secondaryFiles:
  - .tbi
  type:
    items: File
    type: array
- id: config__algorithm__coverage_interval
  type:
    items:
    - 'null'
    - string
    type: array
- id: genome_resources__rnaseq__gene_bed
  type:
    items: File
    type: array
- id: genome_resources__variation__train_hapmap
  secondaryFiles:
  - .tbi
  type:
    items: File
    type: array
- id: rgnames__lb
  type:
    items:
    - 'null'
    - string
    type: array
- id: genome_resources__variation__clinvar
  secondaryFiles:
  - .tbi
  type:
    items: File
    type: array
- id: genome_resources__variation__esp
  secondaryFiles:
  - .tbi
  type:
    items: File
    type: array
- id: rgnames__rg
  type:
    items: string
    type: array
- id: metadata__batch
  type:
    items: string
    type: array
- id: rgnames__lane
  type:
    items: string
    type: array
- id: genome_resources__variation__lcr
  type:
    items:
    - 'null'
    - string
    type: array
- id: config__algorithm__min_allele_fraction
  type:
    items: long
    type: array
- id: config__algorithm__nomap_split_targets
  type:
    items: long
    type: array
- id: reference__bwa__indexes
  secondaryFiles:
  - ^.ann
  - ^.pac
  - ^.sa
  - ^.bwt
  type:
    items: File
    type: array
- id: vrn_file
  type:
    items:
    - 'null'
    - string
    type: array
- id: reference__twobit
  type:
    items: File
    type: array
- id: reference__genome_context
  secondaryFiles:
  - .tbi
  type:
    items:
      items: File
      type: array
    type: array
- id: config__algorithm__bam_clean
  type:
    items:
    - string
    - 'null'
    - boolean
    type: array
- id: config__algorithm__nomap_split_size
  type:
    items: long
    type: array
- id: config__algorithm__validate
  type:
    items:
    - File
    - 'null'
    type: array
- id: reference__snpeff__hg19
  type:
    items: File
    type: array
- id: description
  type:
    items: string
    type: array
- id: config__algorithm__validate_regions
  type:
    items:
    - File
    - 'null'
    type: array
- id: config__algorithm__aligner
  type:
    items: string
    type: array
- id: rgnames__pl
  type:
    items: string
    type: array
- id: genome_build
  type:
    items: string
    type: array
- id: rgnames__pu
  type:
    items: string
    type: array
- id: genome_resources__variation__exac
  secondaryFiles:
  - .tbi
  type:
    items: File
    type: array
- id: config__algorithm__recalibrate
  type:
    items:
    - string
    - 'null'
    - boolean
    type: array
- id: config__algorithm__coverage
  type:
    items:
    - File
    - 'null'
    type: array
- id: genome_resources__aliases__human
  type:
    items:
    - string
    - 'null'
    - boolean
    type: array
- id: config__algorithm__tools_off
  type:
    items:
      items: string
      type: array
    type: array
- id: genome_resources__variation__dbsnp
  secondaryFiles:
  - .tbi
  type:
    items: File
    type: array
- id: config__algorithm__mark_duplicates
  type:
    items:
    - string
    - 'null'
    - boolean
    type: array
- id: genome_resources__variation__polyx
  type:
    items:
    - 'null'
    - string
    type: array
- id: genome_resources__variation__encode_blacklist
  type:
    items:
    - 'null'
    - string
    type: array
- id: genome_resources__variation__cosmic
  secondaryFiles:
  - .tbi
  type:
    items: File
    type: array
- id: config__algorithm__ensemble
  type:
    items:
    - 'null'
    - string
    type: array
- id: config__algorithm__qc
  type:
    items:
      items: string
      type: array
    type: array
- id: analysis
  type:
    items: string
    type: array
- id: rgnames__sample
  type:
    items: string
    type: array
- id: config__algorithm__tools_on
  type:
    items:
      items: string
      type: array
    type: array
- id: config__algorithm__effects
  type:
    items:
    - string
    - 'null'
    - boolean
    type: array
- id: config__algorithm__variant_regions
  type:
    items: File
    type: array
- id: genome_resources__aliases__ensembl
  type:
    items: string
    type: array
- id: config__algorithm__exclude_regions
  type:
    items:
    - 'null'
    - items: 'null'
      type: array
    type: array
- id: reference__rtg
  type:
    items: File
    type: array
- id: genome_resources__variation__train_indels
  secondaryFiles:
  - .tbi
  type:
    items: File
    type: array
- id: genome_resources__aliases__snpeff
  type:
    items: string
    type: array
outputs:
- id: rgnames__sample_out
  outputSource: prep_samples/rgnames__sample
  type:
    items: string
    type: array
- id: align_bam
  outputSource: postprocess_alignment/align_bam
  type:
    items:
    - File
    - 'null'
    type: array
- id: regions__sample_callable
  outputSource: postprocess_alignment/regions__sample_callable
  type:
    items:
    - File
    - 'null'
    type: array
- id: validate__grading_summary
  outputSource: summarize_vc/validate__grading_summary
  type:
    items:
    - File
    - 'null'
    type: array
- id: variants__calls
  outputSource: summarize_vc/variants__calls
  type:
    items:
      items:
      - File
      - 'null'
      type: array
    type: array
- id: variants__gvcf
  outputSource: summarize_vc/variants__gvcf
  type:
    items:
    - 'null'
    - items:
      - File
      - 'null'
      type: array
    type: array
- id: summary__multiqc
  outputSource: multiqc_summary/summary__multiqc
  type:
    items:
    - File
    - 'null'
    type: array
requirements:
- class: EnvVarRequirement
  envDef:
  - envName: MPLCONFIGDIR
    envValue: .
- class: ScatterFeatureRequirement
- class: SubworkflowFeatureRequirement
steps:
- id: alignment_to_rec
  in:
  - id: files
    source: files
  - id: analysis
    source: analysis
  - id: config__algorithm__align_split_size
    source: config__algorithm__align_split_size
  - id: reference__fasta__base
    source: reference__fasta__base
  - id: rgnames__pl
    source: rgnames__pl
  - id: rgnames__sample
    source: rgnames__sample
  - id: rgnames__pu
    source: rgnames__pu
  - id: rgnames__lane
    source: rgnames__lane
  - id: rgnames__rg
    source: rgnames__rg
  - id: rgnames__lb
    source: rgnames__lb
  - id: reference__bwa__indexes
    source: reference__bwa__indexes
  - id: config__algorithm__aligner
    source: config__algorithm__aligner
  - id: config__algorithm__trim_reads
    source: config__algorithm__trim_reads
  - id: config__algorithm__adapters
    source: config__algorithm__adapters
  - id: config__algorithm__bam_clean
    source: config__algorithm__bam_clean
  - id: config__algorithm__mark_duplicates
    source: config__algorithm__mark_duplicates
  - id: resources
    source: resources
  - id: description
    source: description
  out:
  - id: alignment_rec
  run: steps/alignment_to_rec.cwl
- id: alignment
  in:
  - id: alignment_rec
    source: alignment_to_rec/alignment_rec
  out:
  - id: align_bam
  - id: hla__fastq
  - id: work_bam_plus__disc
  - id: work_bam_plus__sr
  run: wf-alignment.cwl
  scatter:
  - alignment_rec
  scatterMethod: dotproduct
- id: prep_samples_to_rec
  in:
  - id: config__algorithm__coverage
    source: config__algorithm__coverage
  - id: rgnames__sample
    source: rgnames__sample
  - id: config__algorithm__variant_regions
    source: config__algorithm__variant_regions
  - id: reference__fasta__base
    source: reference__fasta__base
  - id: resources
    source: resources
  - id: description
    source: description
  out:
  - id: prep_samples_rec
  run: steps/prep_samples_to_rec.cwl
- id: prep_samples
  in:
  - id: prep_samples_rec
    source: prep_samples_to_rec/prep_samples_rec
  out:
  - id: rgnames__sample
  - id: config__algorithm__variant_regions
  - id: config__algorithm__variant_regions_merged
  - id: config__algorithm__variant_regions_orig
  - id: config__algorithm__coverage
  - id: config__algorithm__coverage_merged
  - id: config__algorithm__coverage_orig
  - id: config__algorithm__seq2c_bed_ready
  run: steps/prep_samples.cwl
  scatter:
  - prep_samples_rec
  scatterMethod: dotproduct
- id: postprocess_alignment_to_rec
  in:
  - id: align_bam
    source: alignment/align_bam
  - id: config__algorithm__coverage_interval
    source: config__algorithm__coverage_interval
  - id: config__algorithm__exclude_regions
    source: config__algorithm__exclude_regions
  - id: config__algorithm__variant_regions
    source: prep_samples/config__algorithm__variant_regions
  - id: config__algorithm__variant_regions_merged
    source: prep_samples/config__algorithm__variant_regions_merged
  - id: config__algorithm__variant_regions_orig
    source: prep_samples/config__algorithm__variant_regions_orig
  - id: config__algorithm__coverage
    source: prep_samples/config__algorithm__coverage
  - id: config__algorithm__coverage_merged
    source: prep_samples/config__algorithm__coverage_merged
  - id: config__algorithm__coverage_orig
    source: prep_samples/config__algorithm__coverage_orig
  - id: config__algorithm__seq2c_bed_ready
    source: prep_samples/config__algorithm__seq2c_bed_ready
  - id: config__algorithm__recalibrate
    source: config__algorithm__recalibrate
  - id: config__algorithm__tools_on
    source: config__algorithm__tools_on
  - id: genome_resources__rnaseq__gene_bed
    source: genome_resources__rnaseq__gene_bed
  - id: genome_resources__variation__dbsnp
    source: genome_resources__variation__dbsnp
  - id: genome_resources__variation__lcr
    source: genome_resources__variation__lcr
  - id: genome_resources__variation__polyx
    source: genome_resources__variation__polyx
  - id: genome_resources__variation__encode_blacklist
    source: genome_resources__variation__encode_blacklist
  - id: reference__twobit
    source: reference__twobit
  - id: reference__fasta__base
    source: reference__fasta__base
  - id: resources
    source: resources
  - id: description
    source: description
  out:
  - id: postprocess_alignment_rec
  run: steps/postprocess_alignment_to_rec.cwl
- id: postprocess_alignment
  in:
  - id: postprocess_alignment_rec
    source: postprocess_alignment_to_rec/postprocess_alignment_rec
  out:
  - id: config__algorithm__coverage_interval
  - id: config__algorithm__variant_regions
  - id: config__algorithm__variant_regions_merged
  - id: config__algorithm__variant_regions_orig
  - id: config__algorithm__coverage
  - id: config__algorithm__coverage_merged
  - id: config__algorithm__coverage_orig
  - id: config__algorithm__seq2c_bed_ready
  - id: regions__callable
  - id: regions__sample_callable
  - id: regions__nblock
  - id: depth__samtools__stats
  - id: depth__samtools__idxstats
  - id: depth__variant_regions__regions
  - id: depth__variant_regions__dist
  - id: depth__sv_regions__regions
  - id: depth__sv_regions__dist
  - id: depth__coverage__regions
  - id: depth__coverage__dist
  - id: depth__coverage__thresholds
  - id: align_bam
  run: steps/postprocess_alignment.cwl
  scatter:
  - postprocess_alignment_rec
  scatterMethod: dotproduct
- id: combine_sample_regions
  in:
  - id: regions__callable
    source: postprocess_alignment/regions__callable
  - id: regions__nblock
    source: postprocess_alignment/regions__nblock
  - id: metadata__batch
    source: metadata__batch
  - id: config__algorithm__nomap_split_size
    source: config__algorithm__nomap_split_size
  - id: config__algorithm__nomap_split_targets
    source: config__algorithm__nomap_split_targets
  - id: reference__fasta__base
    source: reference__fasta__base
  - id: resources
    source: resources
  - id: description
    source: description
  out:
  - id: config__algorithm__callable_regions
  - id: config__algorithm__non_callable_regions
  - id: config__algorithm__callable_count
  run: steps/combine_sample_regions.cwl
- id: batch_for_variantcall
  in:
  - id: analysis
    source: analysis
  - id: genome_build
    source: genome_build
  - id: align_bam
    source: postprocess_alignment/align_bam
  - id: vrn_file
    source: vrn_file
  - id: metadata__batch
    source: metadata__batch
  - id: metadata__phenotype
    source: metadata__phenotype
  - id: config__algorithm__callable_regions
    source: combine_sample_regions/config__algorithm__callable_regions
  - id: regions__sample_callable
    source: postprocess_alignment/regions__sample_callable
  - id: config__algorithm__variantcaller
    source: config__algorithm__variantcaller
  - id: config__algorithm__ensemble
    source: config__algorithm__ensemble
  - id: config__algorithm__vcfanno
    source: config__algorithm__vcfanno
  - id: config__algorithm__coverage_interval
    source: postprocess_alignment/config__algorithm__coverage_interval
  - id: config__algorithm__effects
    source: config__algorithm__effects
  - id: config__algorithm__min_allele_fraction
    source: config__algorithm__min_allele_fraction
  - id: config__algorithm__exclude_regions
    source: config__algorithm__exclude_regions
  - id: config__algorithm__variant_regions
    source: postprocess_alignment/config__algorithm__variant_regions
  - id: config__algorithm__variant_regions_merged
    source: postprocess_alignment/config__algorithm__variant_regions_merged
  - id: config__algorithm__validate
    source: config__algorithm__validate
  - id: config__algorithm__validate_regions
    source: config__algorithm__validate_regions
  - id: config__algorithm__tools_on
    source: config__algorithm__tools_on
  - id: config__algorithm__tools_off
    source: config__algorithm__tools_off
  - id: reference__fasta__base
    source: reference__fasta__base
  - id: reference__twobit
    source: reference__twobit
  - id: reference__rtg
    source: reference__rtg
  - id: reference__genome_context
    source: reference__genome_context
  - id: genome_resources__variation__clinvar
    source: genome_resources__variation__clinvar
  - id: genome_resources__variation__cosmic
    source: genome_resources__variation__cosmic
  - id: genome_resources__variation__dbsnp
    source: genome_resources__variation__dbsnp
  - id: genome_resources__variation__esp
    source: genome_resources__variation__esp
  - id: genome_resources__variation__exac
    source: genome_resources__variation__exac
  - id: genome_resources__variation__1000g
    source: genome_resources__variation__1000g
  - id: genome_resources__variation__lcr
    source: genome_resources__variation__lcr
  - id: genome_resources__variation__polyx
    source: genome_resources__variation__polyx
  - id: genome_resources__variation__encode_blacklist
    source: genome_resources__variation__encode_blacklist
  - id: genome_resources__aliases__ensembl
    source: genome_resources__aliases__ensembl
  - id: genome_resources__aliases__human
    source: genome_resources__aliases__human
  - id: genome_resources__aliases__snpeff
    source: genome_resources__aliases__snpeff
  - id: reference__snpeff__hg19
    source: reference__snpeff__hg19
  - id: genome_resources__variation__train_hapmap
    source: genome_resources__variation__train_hapmap
  - id: genome_resources__variation__train_indels
    source: genome_resources__variation__train_indels
  - id: resources
    source: resources
  - id: description
    source: description
  out:
  - id: batch_rec
  run: steps/batch_for_variantcall.cwl
- id: variantcall
  in:
  - id: batch_rec
    source: batch_for_variantcall/batch_rec
  out:
  - id: vc_rec
  run: wf-variantcall.cwl
  scatter:
  - batch_rec
  scatterMethod: dotproduct
- id: batch_for_jointvc
  in:
  - id: vc_rec
    source: variantcall/vc_rec
  out:
  - id: jointvc_batch_rec
  run: steps/batch_for_jointvc.cwl
- id: jointcall
  in:
  - id: jointvc_batch_rec
    source: batch_for_jointvc/jointvc_batch_rec
  out:
  - id: jointvc_rec
  - id: vrn_file_joint
  run: wf-jointcall.cwl
  scatter:
  - jointvc_batch_rec
  scatterMethod: dotproduct
- id: summarize_vc
  in:
  - id: jointvc_rec
    source: jointcall/jointvc_rec
  out:
  - id: variants__calls
  - id: variants__gvcf
  - id: variants__samples
  - id: validate__grading_summary
  - id: validate__grading_plots
  run: steps/summarize_vc.cwl
- id: qc_to_rec
  in:
  - id: align_bam
    source: postprocess_alignment/align_bam
  - id: analysis
    source: analysis
  - id: reference__fasta__base
    source: reference__fasta__base
  - id: config__algorithm__tools_on
    source: config__algorithm__tools_on
  - id: config__algorithm__tools_off
    source: config__algorithm__tools_off
  - id: genome_build
    source: genome_build
  - id: config__algorithm__qc
    source: config__algorithm__qc
  - id: metadata__batch
    source: metadata__batch
  - id: config__algorithm__coverage_interval
    source: postprocess_alignment/config__algorithm__coverage_interval
  - id: depth__variant_regions__regions
    source: postprocess_alignment/depth__variant_regions__regions
  - id: depth__variant_regions__dist
    source: postprocess_alignment/depth__variant_regions__dist
  - id: depth__samtools__stats
    source: postprocess_alignment/depth__samtools__stats
  - id: depth__samtools__idxstats
    source: postprocess_alignment/depth__samtools__idxstats
  - id: depth__sv_regions__regions
    source: postprocess_alignment/depth__sv_regions__regions
  - id: depth__sv_regions__dist
    source: postprocess_alignment/depth__sv_regions__dist
  - id: depth__coverage__regions
    source: postprocess_alignment/depth__coverage__regions
  - id: depth__coverage__dist
    source: postprocess_alignment/depth__coverage__dist
  - id: depth__coverage__thresholds
    source: postprocess_alignment/depth__coverage__thresholds
  - id: config__algorithm__variant_regions
    source: postprocess_alignment/config__algorithm__variant_regions
  - id: config__algorithm__variant_regions_merged
    source: postprocess_alignment/config__algorithm__variant_regions_merged
  - id: config__algorithm__coverage
    source: postprocess_alignment/config__algorithm__coverage
  - id: config__algorithm__coverage_merged
    source: postprocess_alignment/config__algorithm__coverage_merged
  - id: variants__samples
    source: summarize_vc/variants__samples
  - id: resources
    source: resources
  - id: description
    source: description
  out:
  - id: qc_rec
  run: steps/qc_to_rec.cwl
- id: pipeline_summary
  in:
  - id: qc_rec
    source: qc_to_rec/qc_rec
  out:
  - id: qcout_rec
  run: steps/pipeline_summary.cwl
  scatter:
  - qc_rec
  scatterMethod: dotproduct
- id: multiqc_summary
  in:
  - id: qcout_rec
    source: pipeline_summary/qcout_rec
  out:
  - id: summary__multiqc
  run: steps/multiqc_summary.cwl
