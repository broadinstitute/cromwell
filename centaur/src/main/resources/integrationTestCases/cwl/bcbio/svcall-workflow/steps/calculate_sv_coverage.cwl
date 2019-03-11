$namespaces:
  dx: https://www.dnanexus.com/cwl#
arguments:
- position: 0
  valueFrom: sentinel_runtime=cores,$(runtime['cores']),ram,$(runtime['ram'])
- sentinel_parallel=multi-parallel
- sentinel_outputs=sv_rawcoverage_rec:depth__bins__target;depth__bins__antitarget;resources;description;regions__bins__target;regions__bins__antitarget;regions__bins__gcannotated;regions__bins__group;reference__fasta__base;config__algorithm__svcaller;config__algorithm__coverage_interval;genome_resources__rnaseq__gene_bed;metadata__batch;genome_resources__variation__lcr;metadata__phenotype;genome_resources__variation__polyx;genome_resources__variation__encode_blacklist;config__algorithm__sv_regions;config__algorithm__variant_regions;config__algorithm__exclude_regions;align_bam;config__algorithm__variant_regions_merged;depth__variant_regions__regions;config__algorithm__callable_regions
- sentinel_inputs=sv_bin_rec:record
- run_number=0
baseCommand:
- bcbio_nextgen.py
- runfn
- calculate_sv_coverage
- cwl
class: CommandLineTool
cwlVersion: v1.0
hints:
- class: DockerRequirement
  dockerImageId: quay.io/bcbio/bcbio-vc:1.1.4a-741877e
  dockerPull: quay.io/bcbio/bcbio-vc:1.1.4a-741877e
- class: ResourceRequirement
  coresMin: 2
  outdirMin: 1028
  ramMin: 6144
  tmpdirMin: 2
- class: dx:InputResourceRequirement
  indirMin: 1
- class: SoftwareRequirement
  packages:
  - package: mosdepth
    specs:
    - https://anaconda.org/bioconda/mosdepth
  - package: cnvkit
    specs:
    - https://anaconda.org/bioconda/cnvkit
inputs:
- id: sv_bin_rec
  type:
    fields:
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
    - name: resources
      type: string
    - name: description
      type: string
    - name: reference__fasta__base
      type: File
    - name: config__algorithm__svcaller
      type:
        items: string
        type: array
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
    - name: config__algorithm__sv_regions
      type: File
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
    - name: depth__variant_regions__regions
      type:
      - File
      - 'null'
    - name: config__algorithm__callable_regions
      type: File
    name: sv_bin_rec
    type: record
outputs:
- id: sv_rawcoverage_rec
  type:
    fields:
    - name: depth__bins__target
      type:
      - File
      - 'null'
    - name: depth__bins__antitarget
      type:
      - File
      - 'null'
    - name: resources
      type: string
    - name: description
      type: string
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
      type:
        items: string
        type: array
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
    - name: config__algorithm__sv_regions
      type: File
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
    - name: depth__variant_regions__regions
      type:
      - File
      - 'null'
    - name: config__algorithm__callable_regions
      type: File
    name: sv_rawcoverage_rec
    type: record
requirements:
- class: InlineJavascriptRequirement
- class: InitialWorkDirRequirement
  listing:
  - entry: $(JSON.stringify(inputs))
    entryname: cwl.inputs.json
