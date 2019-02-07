class: Workflow
cwlVersion: v1.0
hints: []
inputs:
- id: jointvc_batch_rec
  type:
    items:
      fields:
      - name: resources
        type: string
      - name: description
        type: string
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
      - name: vrn_file
        type: File
      - name: reference__fasta__base
        type: File
      - name: metadata__phenotype
        type: string
      - name: config__algorithm__vcfanno
        type:
          items: File
          type: array
      - name: config__algorithm__variantcaller
        type: string
      - name: genome_resources__variation__1000g
        type: File
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
        type:
        - 'null'
        - string
      - name: config__algorithm__min_allele_fraction
        type: double
      - name: genome_resources__variation__train_hapmap
        type: File
      - name: reference__genome_context
        type:
          items: File
          type: array
      - name: config__algorithm__validate
        type:
        - File
        - 'null'
      - name: reference__snpeff__hg19
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
        type: File
      - name: genome_resources__aliases__human
        type:
        - string
        - 'null'
        - boolean
      - name: config__algorithm__tools_off
        type:
        - 'null'
        - items: 'null'
          type: array
      - name: genome_resources__variation__dbsnp
        type: File
      - name: genome_resources__variation__polyx
        type:
        - 'null'
        - string
      - name: genome_resources__variation__encode_blacklist
        type:
        - 'null'
        - string
      - name: genome_resources__variation__cosmic
        type: File
      - name: config__algorithm__ensemble
        type:
        - 'null'
        - string
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
      - name: reference__rtg
        type: File
      - name: genome_resources__variation__train_indels
        type: File
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
      name: jointvc_batch_rec
      type: record
    type: array
outputs:
- id: jointvc_rec
  outputSource: finalize_jointvc/jointvc_rec
  type:
    items:
      fields:
      - name: resources
        type: string
      - name: description
        type: string
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
      - name: vrn_file
        type: File
      - name: reference__fasta__base
        type: File
      - name: metadata__phenotype
        type: string
      - name: config__algorithm__vcfanno
        type:
          items: File
          type: array
      - name: config__algorithm__variantcaller
        type: string
      - name: genome_resources__variation__1000g
        type: File
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
        type:
        - 'null'
        - string
      - name: config__algorithm__min_allele_fraction
        type: double
      - name: genome_resources__variation__train_hapmap
        type: File
      - name: reference__genome_context
        type:
          items: File
          type: array
      - name: config__algorithm__validate
        type:
        - File
        - 'null'
      - name: reference__snpeff__hg19
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
        type: File
      - name: genome_resources__aliases__human
        type:
        - string
        - 'null'
        - boolean
      - name: config__algorithm__tools_off
        type:
        - 'null'
        - items: 'null'
          type: array
      - name: genome_resources__variation__dbsnp
        type: File
      - name: genome_resources__variation__polyx
        type:
        - 'null'
        - string
      - name: genome_resources__variation__encode_blacklist
        type:
        - 'null'
        - string
      - name: genome_resources__variation__cosmic
        type: File
      - name: config__algorithm__ensemble
        type:
        - 'null'
        - string
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
      - name: reference__rtg
        type: File
      - name: genome_resources__variation__train_indels
        type: File
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
      - name: vrn_file_joint
        type: File
      name: jointvc_rec
      type: record
    type: array
- id: vrn_file_joint
  outputSource: postprocess_variants/vrn_file_joint
  secondaryFiles:
  - .tbi
  type: File
requirements:
- class: EnvVarRequirement
  envDef:
  - envName: MPLCONFIGDIR
    envValue: .
- class: ScatterFeatureRequirement
- class: SubworkflowFeatureRequirement
steps:
- id: get_parallel_regions_jointvc
  in:
  - id: jointvc_batch_rec
    source: jointvc_batch_rec
  out:
  - id: region
  run: steps/get_parallel_regions_jointvc.cwl
- id: run_jointvc
  in:
  - id: jointvc_batch_rec
    source: jointvc_batch_rec
  - id: region_toolinput
    source: get_parallel_regions_jointvc/region
  out:
  - id: vrn_file_region
  - id: region
  run: steps/run_jointvc.cwl
  scatter:
  - region_toolinput
  scatterMethod: dotproduct
- id: concat_batch_variantcalls_jointvc
  in:
  - id: jointvc_batch_rec
    source: jointvc_batch_rec
  - id: region
    source: run_jointvc/region
  - id: vrn_file_region
    source: run_jointvc/vrn_file_region
  out:
  - id: vrn_file_joint
  run: steps/concat_batch_variantcalls_jointvc.cwl
- id: postprocess_variants
  in:
  - id: jointvc_batch_rec
    source: jointvc_batch_rec
  - id: vrn_file_joint_toolinput
    source: concat_batch_variantcalls_jointvc/vrn_file_joint
  out:
  - id: vrn_file_joint
  run: steps/postprocess_variants.cwl
- id: finalize_jointvc
  in:
  - id: jointvc_batch_rec
    source: jointvc_batch_rec
  - id: vrn_file_joint
    source: postprocess_variants/vrn_file_joint
  out:
  - id: jointvc_rec
  run: steps/finalize_jointvc.cwl
