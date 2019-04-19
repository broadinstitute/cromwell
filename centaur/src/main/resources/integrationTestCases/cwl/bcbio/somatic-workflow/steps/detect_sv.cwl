$namespaces:
  arv: http://arvados.org/cwl#
  dx: https://www.dnanexus.com/cwl#
arguments:
- position: 0
  valueFrom: sentinel_runtime=cores,$(runtime['cores']),ram,$(runtime['ram'])
- sentinel_parallel=batch-single
- sentinel_outputs=sv_rec:sv__variantcaller;sv__vrn_file;sv__supplemental;svvalidate__summary;resources;description;config__algorithm__svprioritize;genome_resources__variation__gc_profile;genome_build;genome_resources__variation__germline_het_pon;config__algorithm__tools_off;analysis;config__algorithm__tools_on;config__algorithm__svvalidate;genome_resources__aliases__snpeff;regions__sample_callable;variants__samples;depth__bins__normalized;depth__bins__background;depth__bins__target;depth__bins__antitarget;depth__bins__seq2c;regions__bins__target;regions__bins__antitarget;regions__bins__gcannotated;regions__bins__group;reference__fasta__base;config__algorithm__svcaller;config__algorithm__coverage_interval;genome_resources__rnaseq__gene_bed;metadata__batch;genome_resources__variation__lcr;metadata__phenotype;genome_resources__variation__polyx;genome_resources__variation__encode_blacklist;config__algorithm__variant_regions;config__algorithm__exclude_regions;config__algorithm__variant_regions_merged;config__algorithm__seq2c_bed_ready;depth__variant_regions__regions;config__algorithm__callable_regions
- sentinel_inputs=sv_batch_rec:record
- run_number=0
baseCommand:
- bcbio_nextgen.py
- runfn
- detect_sv
- cwl
class: CommandLineTool
cwlVersion: v1.0
hints:
- class: DockerRequirement
  dockerImageId: quay.io/bcbio/bcbio-vc:1.1.4a-741877e
  dockerPull: quay.io/bcbio/bcbio-vc:1.1.4a-741877e
- class: ResourceRequirement
  coresMin: 2
  outdirMin: 10246
  ramMin: 6144
  tmpdirMin: 3
- class: dx:InputResourceRequirement
  indirMin: 1
- class: SoftwareRequirement
  packages:
  - package: bedtools
    specs:
    - https://anaconda.org/bioconda/bedtools
  - package: cnvkit
    specs:
    - https://anaconda.org/bioconda/cnvkit
  - package: delly
    specs:
    - https://anaconda.org/bioconda/delly
  - package: duphold
    specs:
    - https://anaconda.org/bioconda/duphold
  - package: extract-sv-reads
    specs:
    - https://anaconda.org/bioconda/extract-sv-reads
  - package: lumpy-sv
    specs:
    - https://anaconda.org/bioconda/lumpy-sv
  - package: manta
    specs:
    - https://anaconda.org/bioconda/manta
  - package: break-point-inspector
    specs:
    - https://anaconda.org/bioconda/break-point-inspector
  - package: mosdepth
    specs:
    - https://anaconda.org/bioconda/mosdepth
  - package: samtools
    specs:
    - https://anaconda.org/bioconda/samtools
  - package: smoove
    specs:
    - https://anaconda.org/bioconda/smoove
  - package: pysam>
    specs:
    - https://anaconda.org/bioconda/pysam>
    version:
    - 0.13.0
  - package: seq2c
    specs:
    - https://anaconda.org/bioconda/seq2c
  - package: simple_sv_annotation
    specs:
    - https://anaconda.org/bioconda/simple_sv_annotation
  - package: survivor
    specs:
    - https://anaconda.org/bioconda/survivor
  - package: svtools
    specs:
    - https://anaconda.org/bioconda/svtools
  - package: svtyper
    specs:
    - https://anaconda.org/bioconda/svtyper
  - package: r
    specs:
    - https://anaconda.org/bioconda/r
    version:
    - 3.4.1
  - package: r-base=3.4.1=h4fe35fd_8
    specs:
    - https://anaconda.org/bioconda/r-base=3.4.1=h4fe35fd_8
  - package: xorg-libxt
    specs:
    - https://anaconda.org/bioconda/xorg-libxt
  - package: vawk
    specs:
    - https://anaconda.org/bioconda/vawk
- class: arv:RuntimeConstraints
  keep_cache: 4096
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
- class: InlineJavascriptRequirement
- class: InitialWorkDirRequirement
  listing:
  - entry: $(JSON.stringify(inputs))
    entryname: cwl.inputs.json
