class: Workflow
cwlVersion: v1.0
hints: []
inputs:
- id: sv_batch_rec
  type:
    items:
      fields:
      - name: resources
        type: string
      - name: description
        type: string
      - name: config__algorithm__svprioritize
        type:
        - 'null'
        - string
      - name: genome_resources__variation__gc_profile
        type: File
      - name: reference__snpeff__hg19
        type: File
      - name: genome_build
        type: string
      - name: genome_resources__variation__germline_het_pon
        type:
        - 'null'
        - string
      - name: config__algorithm__tools_off
        type:
        - 'null'
        - items: 'null'
          type: array
      - name: analysis
        type: string
      - name: config__algorithm__tools_on
        type:
        - 'null'
        - items: 'null'
          type: array
      - name: config__algorithm__svvalidate
        type:
        - 'null'
        - string
      - name: genome_resources__aliases__snpeff
        type: string
      - name: work_bam_plus__disc
        type:
        - File
        - 'null'
      - name: work_bam_plus__sr
        type:
        - File
        - 'null'
      - name: regions__sample_callable
        type:
        - File
        - 'null'
      - name: variants__samples
        type:
          items:
            items:
            - File
            - 'null'
            type: array
          type: array
      - name: depth__bins__normalized
        type:
        - File
        - 'null'
      - name: depth__bins__background
        type:
        - File
        - 'null'
      - name: depth__bins__target
        type:
        - File
        - 'null'
      - name: depth__bins__antitarget
        type:
        - File
        - 'null'
      - name: depth__bins__seq2c
        type:
        - File
        - 'null'
      - name: regions__bins__target
        type:
        - File
        - 'null'
      - name: regions__bins__antitarget
        type:
        - File
        - 'null'
      - name: regions__bins__gcannotated
        type:
        - File
        - 'null'
      - name: regions__bins__group
        type:
        - string
        - 'null'
      - name: reference__fasta__base
        type: File
      - name: config__algorithm__svcaller
        type: string
      - name: config__algorithm__coverage_interval
        type:
        - string
        - 'null'
      - name: genome_resources__rnaseq__gene_bed
        type: File
      - name: metadata__batch
        type: string
      - name: genome_resources__variation__lcr
        type:
        - 'null'
        - string
      - name: metadata__phenotype
        type: string
      - name: genome_resources__variation__polyx
        type:
        - 'null'
        - string
      - name: genome_resources__variation__encode_blacklist
        type:
        - 'null'
        - string
      - name: config__algorithm__variant_regions
        type:
        - File
        - 'null'
      - name: config__algorithm__exclude_regions
        type:
        - 'null'
        - items: 'null'
          type: array
      - name: align_bam
        type:
        - File
        - 'null'
      - name: config__algorithm__variant_regions_merged
        type:
        - File
        - 'null'
      - name: config__algorithm__seq2c_bed_ready
        type:
        - File
        - 'null'
      - name: depth__variant_regions__regions
        type:
        - File
        - 'null'
      - name: config__algorithm__callable_regions
        type: File
      name: sv_batch_rec
      type: record
    type: array
outputs:
- id: sv_rec
  outputSource: detect_sv/sv_rec
  type:
    items:
      fields:
      - name: sv__variantcaller
        type:
        - string
        - 'null'
      - name: sv__vrn_file
        type:
        - File
        - 'null'
      - name: sv__supplemental
        type:
          items:
          - File
          type: array
      - name: svvalidate__summary
        type:
        - File
        - 'null'
      - name: resources
        type: string
      - name: description
        type: string
      - name: config__algorithm__svprioritize
        type:
        - 'null'
        - string
      - name: genome_resources__variation__gc_profile
        type: File
      - name: genome_build
        type: string
      - name: genome_resources__variation__germline_het_pon
        type:
        - 'null'
        - string
      - name: config__algorithm__tools_off
        type:
        - 'null'
        - items: 'null'
          type: array
      - name: analysis
        type: string
      - name: config__algorithm__tools_on
        type:
        - 'null'
        - items: 'null'
          type: array
      - name: config__algorithm__svvalidate
        type:
        - 'null'
        - string
      - name: genome_resources__aliases__snpeff
        type: string
      - name: regions__sample_callable
        type:
        - File
        - 'null'
      - name: variants__samples
        type:
          items:
            items:
            - File
            - 'null'
            type: array
          type: array
      - name: depth__bins__normalized
        type:
        - File
        - 'null'
      - name: depth__bins__background
        type:
        - File
        - 'null'
      - name: depth__bins__target
        type:
        - File
        - 'null'
      - name: depth__bins__antitarget
        type:
        - File
        - 'null'
      - name: depth__bins__seq2c
        type:
        - File
        - 'null'
      - name: regions__bins__target
        type:
        - File
        - 'null'
      - name: regions__bins__antitarget
        type:
        - File
        - 'null'
      - name: regions__bins__gcannotated
        type:
        - File
        - 'null'
      - name: regions__bins__group
        type:
        - string
        - 'null'
      - name: reference__fasta__base
        type: File
      - name: config__algorithm__svcaller
        type: string
      - name: config__algorithm__coverage_interval
        type:
        - string
        - 'null'
      - name: genome_resources__rnaseq__gene_bed
        type: File
      - name: metadata__batch
        type: string
      - name: genome_resources__variation__lcr
        type:
        - 'null'
        - string
      - name: metadata__phenotype
        type: string
      - name: genome_resources__variation__polyx
        type:
        - 'null'
        - string
      - name: genome_resources__variation__encode_blacklist
        type:
        - 'null'
        - string
      - name: config__algorithm__variant_regions
        type:
        - File
        - 'null'
      - name: config__algorithm__exclude_regions
        type:
        - 'null'
        - items: 'null'
          type: array
      - name: config__algorithm__variant_regions_merged
        type:
        - File
        - 'null'
      - name: config__algorithm__seq2c_bed_ready
        type:
        - File
        - 'null'
      - name: depth__variant_regions__regions
        type:
        - File
        - 'null'
      - name: config__algorithm__callable_regions
        type: File
      name: sv_rec
      type: record
    type: array
requirements:
- class: EnvVarRequirement
  envDef:
  - envName: MPLCONFIGDIR
    envValue: .
- class: ScatterFeatureRequirement
- class: SubworkflowFeatureRequirement
steps:
- id: detect_sv
  in:
  - id: sv_batch_rec
    source: sv_batch_rec
  out:
  - id: sv_rec
  run: steps/detect_sv.cwl
