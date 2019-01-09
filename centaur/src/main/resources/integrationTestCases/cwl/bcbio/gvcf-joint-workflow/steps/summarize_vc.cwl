$namespaces:
  dx: https://www.dnanexus.com/cwl#
arguments:
- position: 0
  valueFrom: sentinel_runtime=cores,$(runtime['cores']),ram,$(runtime['ram'])
- sentinel_parallel=multi-combined
- sentinel_outputs=variants__calls,variants__gvcf,variants__samples,validate__grading_summary,validate__grading_plots
- sentinel_inputs=jointvc_rec:record
- run_number=0
baseCommand:
- bcbio_nextgen.py
- runfn
- summarize_vc
- cwl
class: CommandLineTool
cwlVersion: v1.0
hints:
- class: DockerRequirement
  dockerImageId: quay.io/bcbio/bcbio-vc
  dockerPull: quay.io/bcbio/bcbio-vc
- class: ResourceRequirement
  coresMin: 1
  outdirMin: 1025
  ramMin: 2048
  tmpdirMin: 1
- class: dx:InputResourceRequirement
  indirMin: 1
inputs:
- id: jointvc_rec
  type:
    items:
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
    type: array
outputs:
- id: variants__calls
  type:
    items:
      items:
      - File
      - 'null'
      type: array
    type: array
- id: variants__gvcf
  type:
    items:
    - 'null'
    - items:
      - File
      - 'null'
      type: array
    type: array
- id: variants__samples
  type:
    items:
      items:
        items:
        - File
        - 'null'
        type: array
      type: array
    type: array
- id: validate__grading_summary
  type:
    items:
    - File
    - 'null'
    type: array
- id: validate__grading_plots
  type:
    items:
      items:
      - File
      - 'null'
      type: array
    type: array
requirements:
- class: InlineJavascriptRequirement
- class: InitialWorkDirRequirement
  listing:
  - entry: $(JSON.stringify(inputs))
    entryname: cwl.inputs.json
