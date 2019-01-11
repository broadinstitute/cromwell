class: Workflow
cwlVersion: v1.0
hints: []
inputs:
- id: batch_rec
  type:
    items:
      fields:
      - name: resources
        type: string
      - name: description
        type: string
      - name: reference__fasta__base
        type: File
      - name: config__algorithm__vcfanno
        type:
        - 'null'
        - items: 'null'
          type: array
      - name: config__algorithm__variantcaller
        type: string
      - name: config__algorithm__coverage_interval
        type:
        - string
        - 'null'
      - name: genome_resources__variation__clinvar
        type: File
      - name: genome_resources__variation__esp
        type: File
      - name: metadata__batch
        type: string
      - name: genome_resources__variation__lcr
        type: File
      - name: config__algorithm__min_allele_fraction
        type: double
      - name: reference__snpeff__GRCh37_75
        type: File
      - name: vrn_file
        type:
        - 'null'
        - string
      - name: genome_resources__variation__train_hapmap
        type: File
      - name: reference__genome_context
        type:
          items: File
          type: array
      - name: reference__rtg
        type: File
      - name: config__algorithm__validate
        type:
        - File
        - 'null'
      - name: genome_resources__variation__1000g
        type: File
      - name: config__algorithm__validate_regions
        type:
        - File
        - 'null'
      - name: genome_build
        type: string
      - name: genome_resources__variation__exac
        type: File
      - name: genome_resources__variation__gnomad_exome
        type:
        - 'null'
        - string
      - name: metadata__phenotype
        type: string
      - name: genome_resources__aliases__human
        type:
        - string
        - 'null'
        - boolean
      - name: config__algorithm__tools_off
        type:
          items: string
          type: array
      - name: genome_resources__variation__dbsnp
        type: File
      - name: genome_resources__variation__polyx
        type: File
      - name: genome_resources__variation__encode_blacklist
        type: File
      - name: genome_resources__variation__cosmic
        type: File
      - name: config__algorithm__ensemble
        type: string
      - name: analysis
        type: string
      - name: config__algorithm__tools_on
        type:
          items: string
          type: array
      - name: config__algorithm__effects
        type: string
      - name: config__algorithm__variant_regions
        type:
        - File
        - 'null'
      - name: genome_resources__aliases__ensembl
        type: string
      - name: config__algorithm__exclude_regions
        type:
        - 'null'
        - items: 'null'
          type: array
      - name: genome_resources__variation__train_indels
        type: File
      - name: genome_resources__aliases__snpeff
        type: string
      - name: align_bam
        type:
        - File
        - 'null'
      - name: config__algorithm__variant_regions_merged
        type:
        - File
        - 'null'
      - name: regions__sample_callable
        type:
        - File
        - 'null'
      - name: config__algorithm__callable_regions
        type: File
      name: batch_rec
      type: record
    type: array
outputs:
- id: vc_rec
  outputSource: compare_to_rm/vc_rec
  type:
    items:
      fields:
      - name: batch_samples
        type:
        - 'null'
        - items: string
          type: array
      - name: validate__summary
        type:
        - File
        - 'null'
      - name: validate__tp
        type:
        - File
        - 'null'
      - name: validate__fp
        type:
        - File
        - 'null'
      - name: validate__fn
        type:
        - File
        - 'null'
      - name: resources
        type: string
      - name: description
        type: string
      - name: vrn_file
        type: File
      - name: reference__fasta__base
        type: File
      - name: config__algorithm__vcfanno
        type:
        - 'null'
        - items: 'null'
          type: array
      - name: config__algorithm__variantcaller
        type: string
      - name: config__algorithm__coverage_interval
        type:
        - string
        - 'null'
      - name: metadata__batch
        type: string
      - name: config__algorithm__min_allele_fraction
        type: double
      - name: reference__snpeff__GRCh37_75
        type: File
      - name: reference__genome_context
        type:
          items: File
          type: array
      - name: reference__rtg
        type: File
      - name: config__algorithm__validate
        type:
        - File
        - 'null'
      - name: config__algorithm__validate_regions
        type:
        - File
        - 'null'
      - name: genome_build
        type: string
      - name: metadata__phenotype
        type: string
      - name: genome_resources__aliases__human
        type:
        - string
        - 'null'
        - boolean
      - name: config__algorithm__tools_off
        type:
          items: string
          type: array
      - name: config__algorithm__ensemble
        type: string
      - name: analysis
        type: string
      - name: config__algorithm__tools_on
        type:
          items: string
          type: array
      - name: config__algorithm__effects
        type: string
      - name: config__algorithm__variant_regions
        type:
        - File
        - 'null'
      - name: genome_resources__aliases__ensembl
        type: string
      - name: config__algorithm__exclude_regions
        type:
        - 'null'
        - items: 'null'
          type: array
      - name: genome_resources__aliases__snpeff
        type: string
      - name: config__algorithm__variant_regions_merged
        type:
        - File
        - 'null'
      - name: regions__sample_callable
        type:
        - File
        - 'null'
      - name: config__algorithm__callable_regions
        type: File
      name: vc_rec
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
- id: get_parallel_regions
  in:
  - id: batch_rec
    source: batch_rec
  out:
  - id: region_block
  run: steps/get_parallel_regions.cwl
- id: variantcall_batch_region
  in:
  - id: batch_rec
    source: batch_rec
  - id: region_block_toolinput
    source: get_parallel_regions/region_block
  out:
  - id: vrn_file_region
  - id: region_block
  run: steps/variantcall_batch_region.cwl
  scatter:
  - region_block_toolinput
  scatterMethod: dotproduct
- id: concat_batch_variantcalls
  in:
  - id: batch_rec
    source: batch_rec
  - id: region_block
    source: variantcall_batch_region/region_block
  - id: vrn_file_region
    source: variantcall_batch_region/vrn_file_region
  out:
  - id: vrn_file
  run: steps/concat_batch_variantcalls.cwl
- id: postprocess_variants
  in:
  - id: batch_rec
    source: batch_rec
  - id: vrn_file_toolinput
    source: concat_batch_variantcalls/vrn_file
  out:
  - id: vrn_file
  run: steps/postprocess_variants.cwl
- id: compare_to_rm
  in:
  - id: batch_rec
    source: batch_rec
  - id: vrn_file
    source: postprocess_variants/vrn_file
  out:
  - id: vc_rec
  run: steps/compare_to_rm.cwl
