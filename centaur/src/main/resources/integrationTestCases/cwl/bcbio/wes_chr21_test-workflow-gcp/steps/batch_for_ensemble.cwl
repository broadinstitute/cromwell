$namespaces:
  dx: https://www.dnanexus.com/cwl#
arguments:
- position: 0
  valueFrom: sentinel_runtime=cores,$(runtime['cores']),ram,$(runtime['ram'])
- sentinel_parallel=multi-combined
- sentinel_outputs=ensemble_prep_rec:batch_id;variants__calls;variants__variantcallers;resources;description;batch_samples;validate__summary;validate__tp;validate__fp;validate__fn;vrn_file;reference__fasta__base;config__algorithm__vcfanno;config__algorithm__variantcaller;config__algorithm__coverage_interval;metadata__batch;config__algorithm__min_allele_fraction;reference__snpeff__GRCh37_75;reference__genome_context;reference__rtg;config__algorithm__validate;config__algorithm__validate_regions;genome_build;metadata__phenotype;genome_resources__aliases__human;config__algorithm__tools_off;config__algorithm__ensemble;analysis;config__algorithm__tools_on;config__algorithm__effects;config__algorithm__variant_regions;genome_resources__aliases__ensembl;config__algorithm__exclude_regions;genome_resources__aliases__snpeff;config__algorithm__variant_regions_merged;regions__sample_callable;config__algorithm__callable_regions
- sentinel_inputs=vc_rec:record
- run_number=0
baseCommand:
- bcbio_nextgen.py
- runfn
- batch_for_ensemble
- cwl
class: CommandLineTool
cwlVersion: v1.0
hints:
- class: DockerRequirement
  dockerImageId: quay.io/bcbio/bcbio-vc:1.1.4a-741877e
  dockerPull: quay.io/bcbio/bcbio-vc:1.1.4a-741877e
- class: ResourceRequirement
  coresMin: 1
  outdirMin: 10240
  ramMin: 3072
  tmpdirMin: 0
- class: dx:InputResourceRequirement
  indirMin: 0
inputs:
- id: vc_rec
  type:
    items:
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
    type: array
outputs:
- id: ensemble_prep_rec
  type:
    items:
      fields:
      - name: batch_id
        type: string
      - name: variants__calls
        type:
          items: File
          type: array
      - name: variants__variantcallers
        type:
          items: string
          type: array
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
      name: ensemble_prep_rec
      type: record
    type: array
requirements:
- class: InlineJavascriptRequirement
- class: InitialWorkDirRequirement
  listing:
  - entry: $(JSON.stringify(inputs))
    entryname: cwl.inputs.json
