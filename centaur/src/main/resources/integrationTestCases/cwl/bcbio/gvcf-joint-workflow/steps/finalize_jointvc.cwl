$namespaces:
  dx: https://www.dnanexus.com/cwl#
arguments:
- position: 0
  valueFrom: sentinel_runtime=cores,$(runtime['cores']),ram,$(runtime['ram'])
- sentinel_parallel=batch-single
- sentinel_outputs=jointvc_rec:resources;description;batch_samples;validate__summary;validate__tp;validate__fp;validate__fn;vrn_file;reference__fasta__base;metadata__phenotype;config__algorithm__vcfanno;config__algorithm__variantcaller;genome_resources__variation__1000g;config__algorithm__coverage_interval;genome_resources__variation__clinvar;genome_resources__variation__esp;metadata__batch;genome_resources__variation__lcr;config__algorithm__min_allele_fraction;genome_resources__variation__train_hapmap;reference__genome_context;config__algorithm__validate;reference__snpeff__hg19;config__algorithm__validate_regions;genome_build;genome_resources__variation__exac;genome_resources__variation__gnomad_exome;genome_resources__aliases__human;config__algorithm__tools_off;genome_resources__variation__dbsnp;genome_resources__variation__polyx;genome_resources__variation__encode_blacklist;genome_resources__variation__cosmic;config__algorithm__ensemble;analysis;config__algorithm__tools_on;config__algorithm__effects;config__algorithm__variant_regions;genome_resources__aliases__ensembl;config__algorithm__exclude_regions;reference__rtg;genome_resources__variation__train_indels;genome_resources__aliases__snpeff;config__algorithm__variant_regions_merged;regions__sample_callable;config__algorithm__callable_regions;vrn_file_joint
- sentinel_inputs=jointvc_batch_rec:record,vrn_file_joint:var
- run_number=0
baseCommand:
- bcbio_nextgen.py
- runfn
- finalize_jointvc
- cwl
class: CommandLineTool
cwlVersion: v1.0
hints:
- class: DockerRequirement
  dockerImageId: quay.io/bcbio/bcbio-vc:1.1.4a-741877e
  dockerPull: quay.io/bcbio/bcbio-vc:1.1.4a-741877e
- class: ResourceRequirement
  coresMin: 1
  outdirMin: 1028
  ramMin: 2048
  tmpdirMin: 2
- class: dx:InputResourceRequirement
  indirMin: 1
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
- id: vrn_file_joint
  secondaryFiles:
  - .tbi
  type: File
outputs:
- id: jointvc_rec
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
requirements:
- class: InlineJavascriptRequirement
- class: InitialWorkDirRequirement
  listing:
  - entry: $(JSON.stringify(inputs))
    entryname: cwl.inputs.json
