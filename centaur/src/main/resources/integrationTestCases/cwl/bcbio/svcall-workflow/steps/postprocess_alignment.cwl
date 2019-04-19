$namespaces:
  arv: http://arvados.org/cwl#
  dx: https://www.dnanexus.com/cwl#
arguments:
- position: 0
  valueFrom: sentinel_runtime=cores,$(runtime['cores']),ram,$(runtime['ram'])
- sentinel_parallel=multi-parallel
- sentinel_outputs=config__algorithm__coverage_interval,config__algorithm__variant_regions,config__algorithm__variant_regions_merged,config__algorithm__variant_regions_orig,config__algorithm__coverage,config__algorithm__coverage_merged,config__algorithm__coverage_orig,config__algorithm__seq2c_bed_ready,regions__callable,regions__sample_callable,regions__nblock,depth__samtools__stats,depth__samtools__idxstats,depth__variant_regions__regions,depth__variant_regions__dist,depth__sv_regions__regions,depth__sv_regions__dist,depth__coverage__regions,depth__coverage__dist,depth__coverage__thresholds,align_bam
- sentinel_inputs=postprocess_alignment_rec:record
- run_number=0
baseCommand:
- bcbio_nextgen.py
- runfn
- postprocess_alignment
- cwl
class: CommandLineTool
cwlVersion: v1.0
hints:
- class: DockerRequirement
  dockerImageId: quay.io/bcbio/bcbio-vc:1.1.4a-741877e
  dockerPull: quay.io/bcbio/bcbio-vc:1.1.4a-741877e
- class: ResourceRequirement
  coresMin: 2
  outdirMin: 1033
  ramMin: 6144
  tmpdirMin: 5
- class: dx:InputResourceRequirement
  indirMin: 1
- class: SoftwareRequirement
  packages:
  - package: sambamba
    specs:
    - https://anaconda.org/bioconda/sambamba
  - package: goleft
    specs:
    - https://anaconda.org/bioconda/goleft
  - package: bedtools
    specs:
    - https://anaconda.org/bioconda/bedtools
  - package: htslib
    specs:
    - https://anaconda.org/bioconda/htslib
  - package: gatk4
    specs:
    - https://anaconda.org/bioconda/gatk4
  - package: mosdepth
    specs:
    - https://anaconda.org/bioconda/mosdepth
  - package: sentieon
    specs:
    - https://anaconda.org/bioconda/sentieon
- class: arv:APIRequirement
inputs:
- id: postprocess_alignment_rec
  type:
    fields:
    - name: resources
      type: string
    - name: description
      type: string
    - name: reference__fasta__base
      type: File
    - name: config__algorithm__coverage_interval
      type:
      - 'null'
      - string
    - name: genome_resources__rnaseq__gene_bed
      type: File
    - name: genome_resources__variation__lcr
      type:
      - 'null'
      - string
    - name: config__algorithm__recalibrate
      type:
      - string
      - 'null'
      - boolean
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
    - name: config__algorithm__tools_on
      type:
        items: string
        type: array
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
    - name: config__algorithm__variant_regions_orig
      type:
      - File
      - 'null'
    - name: config__algorithm__coverage
      type:
      - File
      - 'null'
    - name: config__algorithm__coverage_merged
      type:
      - File
      - 'null'
    - name: config__algorithm__coverage_orig
      type:
      - File
      - 'null'
    - name: config__algorithm__seq2c_bed_ready
      type:
      - File
      - 'null'
    name: postprocess_alignment_rec
    type: record
outputs:
- id: config__algorithm__coverage_interval
  type:
  - string
  - 'null'
- id: config__algorithm__variant_regions
  type:
  - File
  - 'null'
- id: config__algorithm__variant_regions_merged
  type:
  - File
  - 'null'
- id: config__algorithm__variant_regions_orig
  type:
  - File
  - 'null'
- id: config__algorithm__coverage
  type:
  - File
  - 'null'
- id: config__algorithm__coverage_merged
  type:
  - File
  - 'null'
- id: config__algorithm__coverage_orig
  type:
  - File
  - 'null'
- id: config__algorithm__seq2c_bed_ready
  type:
  - File
  - 'null'
- id: regions__callable
  type:
  - File
  - 'null'
- id: regions__sample_callable
  type:
  - File
  - 'null'
- id: regions__nblock
  type:
  - File
  - 'null'
- id: depth__samtools__stats
  type:
  - File
  - 'null'
- id: depth__samtools__idxstats
  type:
  - File
  - 'null'
- id: depth__variant_regions__regions
  type:
  - File
  - 'null'
- id: depth__variant_regions__dist
  type:
  - File
  - 'null'
- id: depth__sv_regions__regions
  type:
  - File
  - 'null'
- id: depth__sv_regions__dist
  type:
  - File
  - 'null'
- id: depth__coverage__regions
  type:
  - File
  - 'null'
- id: depth__coverage__dist
  type:
  - File
  - 'null'
- id: depth__coverage__thresholds
  type:
  - File
  - 'null'
- id: align_bam
  secondaryFiles:
  - .bai
  type:
  - File
  - 'null'
requirements:
- class: InlineJavascriptRequirement
- class: InitialWorkDirRequirement
  listing:
  - entry: $(JSON.stringify(inputs))
    entryname: cwl.inputs.json
