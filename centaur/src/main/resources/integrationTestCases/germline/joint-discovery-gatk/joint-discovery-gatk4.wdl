## Copyright Broad Institute, 2017
##
## This WDL implements the joint discovery and VQSR filtering portion of the GATK
## Best Practices (June 2016) for germline SNP and Indel discovery in human
## whole-genome sequencing (WGS) and exome sequencing data.
##
## Requirements/expectations :
## - One or more GVCFs produced by HaplotypeCaller in GVCF mode
## - Bare minimum 1 WGS sample or 30 Exome samples. Gene panels are not supported.
##
## Outputs :
## - A VCF file and its index, filtered using variant quality score recalibration
##   (VQSR) with genotypes for all samples present in the input VCF. All sites that
##   are present in the input VCF are retained; filtered sites are annotated as such
##   in the FILTER field.
##
## Note about VQSR wiring :
## The SNP and INDEL models are built in parallel, but then the corresponding
## recalibrations are applied in series. Because the INDEL model is generally ready
## first (because there are fewer indels than SNPs) we set INDEL recalibration to
## be applied first to the input VCF, while the SNP model is still being built. By
## the time the SNP model is available, the indel-recalibrated file is available to
## serve as input to apply the SNP recalibration. If we did it the other way around,
## we would have to wait until the SNP recal file was available despite the INDEL
## recal file being there already, then apply SNP recalibration, then apply INDEL
## recalibration. This would lead to a longer wall clock time for complete workflow
## execution. Wiring the INDEL recalibration to be applied first solves the problem.
##
## Cromwell version support
## - Successfully tested on v29
## - Does not work on versions < v23 due to output syntax
##
## Runtime parameters are optimized for Broad's Google Cloud Platform implementation.
## For program versions, see docker containers.
##
## LICENSING :
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.

workflow JointGenotyping {
  File unpadded_intervals_file

  String callset_name
  File sample_name_map

  File ref_fasta
  File ref_fasta_index
  File ref_dict

  File dbsnp_vcf
  File dbsnp_vcf_index

  String gatk_docker
  String gatk_path
  String python_docker

  Int small_disk
  Int medium_disk
  Int huge_disk

  Array[String] snp_recalibration_tranche_values
  Array[String] snp_recalibration_annotation_values
  Array[String] indel_recalibration_tranche_values
  Array[String] indel_recalibration_annotation_values

  File eval_interval_list
  File hapmap_resource_vcf
  File hapmap_resource_vcf_index
  File omni_resource_vcf
  File omni_resource_vcf_index
  File one_thousand_genomes_resource_vcf
  File one_thousand_genomes_resource_vcf_index
  File mills_resource_vcf
  File mills_resource_vcf_index
  File axiomPoly_resource_vcf
  File axiomPoly_resource_vcf_index
  File dbsnp_resource_vcf = dbsnp_vcf
  File dbsnp_resource_vcf_index = dbsnp_vcf_index

  # ExcessHet is a phred-scaled p-value. We want a cutoff of anything more extreme
  # than a z-score of -4.5 which is a p-value of 3.4e-06, which phred-scaled is 54.69
  Float excess_het_threshold = 54.69
  Float snp_filter_level
  Float indel_filter_level
  Int SNP_VQSR_downsampleFactor

  Int num_of_original_intervals = length(read_lines(unpadded_intervals_file))
  Int num_gvcfs = length(read_lines(sample_name_map))

  # Make a 2.5:1 interval number to samples in callset ratio interval list
  Int possible_merge_count = floor(num_of_original_intervals / num_gvcfs / 2.5)
  Int merge_count = if possible_merge_count > 1 then possible_merge_count else 1

  call DynamicallyCombineIntervals {
    input:
      intervals = unpadded_intervals_file,
      merge_count = merge_count,
      docker_image = python_docker
  }

  Array[String] unpadded_intervals = read_lines(DynamicallyCombineIntervals.output_intervals)

  scatter (idx in range(length(unpadded_intervals))) {
    # the batch_size value was carefully chosen here as it
    # is the optimal value for the amount of memory allocated
    # within the task; please do not change it without consulting
    # the Hellbender (GATK engine) team!
    call ImportGVCFs {
      input:
        sample_name_map = sample_name_map,
        interval = unpadded_intervals[idx],
        workspace_dir_name = "genomicsdb",
        disk_size = medium_disk,
        docker_image = gatk_docker,
        gatk_path = gatk_path,
        batch_size = 50
    }

    call GenotypeGVCFs {
      input:
        workspace_tar = ImportGVCFs.output_genomicsdb,
        interval = unpadded_intervals[idx],
        output_vcf_filename = "output.vcf.gz",
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        dbsnp_vcf = dbsnp_vcf,
        disk_size = medium_disk,
        docker_image = gatk_docker,
        gatk_path = gatk_path
    }

    call HardFilterAndMakeSitesOnlyVcf {
      input:
        vcf = GenotypeGVCFs.output_vcf,
        vcf_index = GenotypeGVCFs.output_vcf_index,
        excess_het_threshold = excess_het_threshold,
        variant_filtered_vcf_filename = callset_name + "." + idx + ".variant_filtered.vcf.gz",
        sites_only_vcf_filename = callset_name + "." + idx + ".sites_only.variant_filtered.vcf.gz",
        disk_size = medium_disk,
        docker_image = gatk_docker,
        gatk_path = gatk_path
    }
  }

  call GatherVcfs as SitesOnlyGatherVcf {
    input:
      input_vcfs_fofn = write_lines(HardFilterAndMakeSitesOnlyVcf.sites_only_vcf),
      output_vcf_name = callset_name + ".sites_only.vcf.gz",
      disk_size = medium_disk,
      docker_image = gatk_docker,
      gatk_path = gatk_path
  }

  call SNPsVariantRecalibratorCreateModel {
    input:
      sites_only_variant_filtered_vcf = SitesOnlyGatherVcf.output_vcf,
      sites_only_variant_filtered_vcf_index = SitesOnlyGatherVcf.output_vcf_index,
      recalibration_filename = callset_name + ".snps.recal",
      tranches_filename = callset_name + ".snps.tranches",
      recalibration_tranche_values = snp_recalibration_tranche_values,
      recalibration_annotation_values = snp_recalibration_annotation_values,
      downsampleFactor = SNP_VQSR_downsampleFactor,
      model_report_filename = callset_name + ".snps.model.report",
      hapmap_resource_vcf = hapmap_resource_vcf,
      hapmap_resource_vcf_index = hapmap_resource_vcf_index,
      omni_resource_vcf = omni_resource_vcf,
      omni_resource_vcf_index = omni_resource_vcf_index,
      one_thousand_genomes_resource_vcf = one_thousand_genomes_resource_vcf,
      one_thousand_genomes_resource_vcf_index = one_thousand_genomes_resource_vcf_index,
      dbsnp_resource_vcf = dbsnp_resource_vcf,
      dbsnp_resource_vcf_index = dbsnp_resource_vcf_index,
      disk_size = small_disk,
      docker_image = gatk_docker,
      gatk_path = gatk_path
  }

  call IndelsVariantRecalibrator {
    input:
      sites_only_variant_filtered_vcf = SitesOnlyGatherVcf.output_vcf,
      sites_only_variant_filtered_vcf_index = SitesOnlyGatherVcf.output_vcf_index,
      recalibration_filename = callset_name + ".indels.recal",
      tranches_filename = callset_name + ".indels.tranches",
      recalibration_tranche_values = indel_recalibration_tranche_values,
      recalibration_annotation_values = indel_recalibration_annotation_values,
      mills_resource_vcf = mills_resource_vcf,
      mills_resource_vcf_index = mills_resource_vcf_index,
      axiomPoly_resource_vcf = axiomPoly_resource_vcf,
      axiomPoly_resource_vcf_index = axiomPoly_resource_vcf_index,
      dbsnp_resource_vcf = dbsnp_resource_vcf,
      dbsnp_resource_vcf_index = dbsnp_resource_vcf_index,
      disk_size = small_disk,
      docker_image = gatk_docker,
      gatk_path = gatk_path
  }

  scatter (idx in range(length(HardFilterAndMakeSitesOnlyVcf.sites_only_vcf))) {
    call SNPsVariantRecalibratorScattered {
      input:
        sites_only_variant_filtered_vcf = HardFilterAndMakeSitesOnlyVcf.sites_only_vcf[idx],
        sites_only_variant_filtered_vcf_index = HardFilterAndMakeSitesOnlyVcf.sites_only_vcf_index[idx],
        recalibration_filename = callset_name + ".snps." + idx + ".recal",
        tranches_filename = callset_name + ".snps." + idx + ".tranches",
        recalibration_tranche_values = snp_recalibration_tranche_values,
        recalibration_annotation_values = snp_recalibration_annotation_values,
        model_report = SNPsVariantRecalibratorCreateModel.model_report,
        hapmap_resource_vcf = hapmap_resource_vcf,
        hapmap_resource_vcf_index = hapmap_resource_vcf_index,
        omni_resource_vcf = omni_resource_vcf,
        omni_resource_vcf_index = omni_resource_vcf_index,
        one_thousand_genomes_resource_vcf = one_thousand_genomes_resource_vcf,
        one_thousand_genomes_resource_vcf_index = one_thousand_genomes_resource_vcf_index,
        dbsnp_resource_vcf = dbsnp_resource_vcf,
        dbsnp_resource_vcf_index = dbsnp_resource_vcf_index,
        disk_size = small_disk,
        docker_image = gatk_docker,
        gatk_path = gatk_path
      }
  }

  call GatherTranches as SNPGatherTranches {
    input:
      input_fofn = write_lines(SNPsVariantRecalibratorScattered.tranches),
      output_filename = callset_name + ".snps.gathered.tranches",
      disk_size = small_disk,
      docker_image = gatk_docker,
      gatk_path = gatk_path

  }

  # For small callsets (fewer than 1000 samples) we can gather the VCF shards and collect metrics directly.
  # For anything larger, we need to keep the VCF sharded and gather metrics collected from them.
  Boolean is_small_callset = num_gvcfs <= 1000

  scatter (idx in range(length(HardFilterAndMakeSitesOnlyVcf.variant_filtered_vcf))) {
    call ApplyRecalibration {
      input:
        recalibrated_vcf_filename = callset_name + ".filtered." + idx + ".vcf.gz",
        input_vcf = HardFilterAndMakeSitesOnlyVcf.variant_filtered_vcf[idx],
        input_vcf_index = HardFilterAndMakeSitesOnlyVcf.variant_filtered_vcf_index[idx],
        indels_recalibration = IndelsVariantRecalibrator.recalibration,
        indels_recalibration_index = IndelsVariantRecalibrator.recalibration_index,
        indels_tranches = IndelsVariantRecalibrator.tranches,
        snps_recalibration = SNPsVariantRecalibratorScattered.recalibration[idx],
        snps_recalibration_index = SNPsVariantRecalibratorScattered.recalibration_index[idx],
        snps_tranches = SNPGatherTranches.tranches,
        indel_filter_level = indel_filter_level,
        snp_filter_level = snp_filter_level,
        disk_size = medium_disk,
        docker_image = gatk_docker,
        gatk_path = gatk_path
    }

    # for large callsets we need to collect metrics from the shards and gather them later
    if (!is_small_callset) {
      call CollectVariantCallingMetrics as CollectMetricsSharded {
        input:
          input_vcf = ApplyRecalibration.recalibrated_vcf,
          input_vcf_index = ApplyRecalibration.recalibrated_vcf_index,
          metrics_filename_prefix = callset_name + "." + idx,
          dbsnp_vcf = dbsnp_vcf,
          dbsnp_vcf_index = dbsnp_vcf_index,
          interval_list = eval_interval_list,
          ref_dict = ref_dict,
          disk_size = small_disk,
          docker_image = gatk_docker,
          gatk_path = gatk_path
      }
    }
  }

  # for small callsets we can gather the VCF shards and then collect metrics on it
  if (is_small_callset) {
    call GatherVcfs as FinalGatherVcf {
      input:
        input_vcfs_fofn = write_lines(ApplyRecalibration.recalibrated_vcf),
        output_vcf_name = callset_name + ".vcf.gz",
        disk_size = huge_disk,
        docker_image = gatk_docker,
        gatk_path = gatk_path
    }

    call CollectVariantCallingMetrics as CollectMetricsOnFullVcf {
      input:
        input_vcf = FinalGatherVcf.output_vcf,
        input_vcf_index = FinalGatherVcf.output_vcf_index,
        metrics_filename_prefix = callset_name,
        dbsnp_vcf = dbsnp_vcf,
        dbsnp_vcf_index = dbsnp_vcf_index,
        interval_list = eval_interval_list,
        ref_dict = ref_dict,
        disk_size = small_disk,
        docker_image = gatk_docker,
        gatk_path = gatk_path
    }
  }

  # for large callsets we still need to gather the sharded metrics
  if (!is_small_callset) {
    call GatherMetrics {
      input:
        input_details_fofn = write_lines(select_all(CollectMetricsSharded.detail_metrics_file)),
        input_summaries_fofn = write_lines(select_all(CollectMetricsSharded.summary_metrics_file)),
        output_prefix = callset_name,
        disk_size = medium_disk,
        docker_image = gatk_docker,
        gatk_path = gatk_path
    }
  }

  output {
    # outputs from the small callset path through the wdl
    FinalGatherVcf.output_vcf
    FinalGatherVcf.output_vcf_index
    CollectMetricsOnFullVcf.detail_metrics_file
    CollectMetricsOnFullVcf.summary_metrics_file

    # outputs from the large callset path through the wdl
    # (note that we do not list ApplyRecalibration here because it is run in both paths)
    GatherMetrics.detail_metrics_file
    GatherMetrics.summary_metrics_file

    # output the interval list generated/used by this run workflow
    DynamicallyCombineIntervals.output_intervals
  }
}

task ImportGVCFs {
  File sample_name_map
  String interval

  String workspace_dir_name

  String java_opt
  String gatk_path

  String docker_image
  Int disk_size
  String mem_size
  Int preemptibles
  Int batch_size

  command <<<
    set -e

    rm -rf ${workspace_dir_name}

    # The memory setting here is very important and must be several GB lower
    # than the total memory allocated to the VM because this tool uses
    # a significant amount of non-heap memory for native libraries.
    # Also, testing has shown that the multithreaded reader initialization
    # does not scale well beyond 5 threads, so don't increase beyond that.
    ${gatk_path} --java-options "${java_opt}" \
    GenomicsDBImport \
    --genomicsdb-workspace-path ${workspace_dir_name} \
    --batch-size ${batch_size} \
    -L ${interval} \
    --sample-name-map ${sample_name_map} \
    --reader-threads 5 \
    -ip 500

    tar -cf ${workspace_dir_name}.tar ${workspace_dir_name}

  >>>
  runtime {
    docker: docker_image
    memory: mem_size
    cpu: "2"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: preemptibles
  }
  output {
    File output_genomicsdb = "${workspace_dir_name}.tar"
  }
}

task GenotypeGVCFs {
  File workspace_tar
  String interval

  String output_vcf_filename

  String gatk_path
  String java_opt

  File ref_fasta
  File ref_fasta_index
  File ref_dict

  String dbsnp_vcf

  String docker_image
  Int disk_size
  String mem_size
  Int preemptibles

  command <<<
    set -e

    tar -xf ${workspace_tar}
    WORKSPACE=$( basename ${workspace_tar} .tar)

    ${gatk_path} --java-options "${java_opt}" \
     GenotypeGVCFs \
     -R ${ref_fasta} \
     -O ${output_vcf_filename} \
     -D ${dbsnp_vcf} \
     -G StandardAnnotation \
     --only-output-calls-starting-in-intervals \
     -new-qual \
     -V gendb://$WORKSPACE \
     -L ${interval}
  >>>
  runtime {
    docker: docker_image
    memory: mem_size
    cpu: "2"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: preemptibles
  }
  output {
    File output_vcf = "${output_vcf_filename}"
    File output_vcf_index = "${output_vcf_filename}.tbi"
  }
}

task HardFilterAndMakeSitesOnlyVcf {
  File vcf
  File vcf_index
  Float excess_het_threshold

  String variant_filtered_vcf_filename
  String sites_only_vcf_filename

  String gatk_path
  String java_opt

  String docker_image
  Int disk_size
  String mem_size
  Int preemptibles

  command {
    set -e

    ${gatk_path} --java-options "${java_opt}" \
      VariantFiltration \
      --filter-expression "ExcessHet > ${excess_het_threshold}" \
      --filter-name ExcessHet \
      -O ${variant_filtered_vcf_filename} \
      -V ${vcf}

    ${gatk_path} --java-options "${java_opt}" \
      MakeSitesOnlyVcf \
      --INPUT ${variant_filtered_vcf_filename} \
      --OUTPUT ${sites_only_vcf_filename}

  }
  runtime {
    docker: docker_image
    memory: mem_size
    cpu: "1"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: preemptibles
  }
  output {
    File variant_filtered_vcf = "${variant_filtered_vcf_filename}"
    File variant_filtered_vcf_index = "${variant_filtered_vcf_filename}.tbi"
    File sites_only_vcf = "${sites_only_vcf_filename}"
    File sites_only_vcf_index = "${sites_only_vcf_filename}.tbi"
  }
}

task IndelsVariantRecalibrator {
  String recalibration_filename
  String tranches_filename

  Array[String] recalibration_tranche_values
  Array[String] recalibration_annotation_values

  File sites_only_variant_filtered_vcf
  File sites_only_variant_filtered_vcf_index

  File mills_resource_vcf
  File axiomPoly_resource_vcf
  File dbsnp_resource_vcf
  File mills_resource_vcf_index
  File axiomPoly_resource_vcf_index
  File dbsnp_resource_vcf_index

  String gatk_path
  String java_opt

  String docker_image
  Int disk_size
  String mem_size
  Int preemptibles

  command {
    ${gatk_path} --java-options "${java_opt}" \
      VariantRecalibrator \
      -V ${sites_only_variant_filtered_vcf} \
      -O ${recalibration_filename} \
      --tranches-file ${tranches_filename} \
      -trust-all-polymorphic \
      -tranche ${sep=' -tranche ' recalibration_tranche_values} \
      -an ${sep=' -an ' recalibration_annotation_values} \
      -mode INDEL \
      --max-gaussians 4 \
      -resource mills,known=false,training=true,truth=true,prior=12:${mills_resource_vcf} \
      -resource axiomPoly,known=false,training=true,truth=false,prior=10:${axiomPoly_resource_vcf} \
      -resource dbsnp,known=true,training=false,truth=false,prior=2:${dbsnp_resource_vcf}
  }
  runtime {
    docker: docker_image
    memory: mem_size
    cpu: "2"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: preemptibles
  }
  output {
    File recalibration = "${recalibration_filename}"
    File recalibration_index = "${recalibration_filename}.idx"
    File tranches = "${tranches_filename}"
  }
}

task SNPsVariantRecalibratorCreateModel {
  String recalibration_filename
  String tranches_filename
  Int downsampleFactor
  String model_report_filename

  Array[String] recalibration_tranche_values
  Array[String] recalibration_annotation_values

  File sites_only_variant_filtered_vcf
  File sites_only_variant_filtered_vcf_index

  File hapmap_resource_vcf
  File omni_resource_vcf
  File one_thousand_genomes_resource_vcf
  File dbsnp_resource_vcf
  File hapmap_resource_vcf_index
  File omni_resource_vcf_index
  File one_thousand_genomes_resource_vcf_index
  File dbsnp_resource_vcf_index

  String gatk_path
  String java_opt

  String docker_image
  Int disk_size
  String mem_size
  Int preemptibles

  command {
    ${gatk_path} --java-options "${java_opt}" \
      VariantRecalibrator \
      -V ${sites_only_variant_filtered_vcf} \
      -O ${recalibration_filename} \
      --tranches-file ${tranches_filename} \
      -trust-all-polymorphic \
      -tranche ${sep=' -tranche ' recalibration_tranche_values} \
      -an ${sep=' -an ' recalibration_annotation_values} \
      -mode SNP \
      -sample-every ${downsampleFactor} \
      --output-model ${model_report_filename} \
      --max-gaussians 6 \
      -resource hapmap,known=false,training=true,truth=true,prior=15:${hapmap_resource_vcf} \
      -resource omni,known=false,training=true,truth=true,prior=12:${omni_resource_vcf} \
      -resource 1000G,known=false,training=true,truth=false,prior=10:${one_thousand_genomes_resource_vcf} \
      -resource dbsnp,known=true,training=false,truth=false,prior=7:${dbsnp_resource_vcf}
  }
  runtime {
    docker: docker_image
    memory: mem_size
    cpu: "2"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: preemptibles
  }
  output {
    File model_report = "${model_report_filename}"
  }
}

task SNPsVariantRecalibratorScattered {
  String recalibration_filename
  String tranches_filename
  File model_report

  Array[String] recalibration_tranche_values
  Array[String] recalibration_annotation_values

  File sites_only_variant_filtered_vcf
  File sites_only_variant_filtered_vcf_index

  File hapmap_resource_vcf
  File omni_resource_vcf
  File one_thousand_genomes_resource_vcf
  File dbsnp_resource_vcf
  File hapmap_resource_vcf_index
  File omni_resource_vcf_index
  File one_thousand_genomes_resource_vcf_index
  File dbsnp_resource_vcf_index

  String gatk_path
  String java_opt

  String docker_image
  Int disk_size
  String mem_size
  Int preemptibles

  command {
    ${gatk_path} --java-options "${java_opt}" \
      VariantRecalibrator \
      -V ${sites_only_variant_filtered_vcf} \
      -O ${recalibration_filename} \
      --tranches-file ${tranches_filename} \
      -trust-all-polymorphic \
      -tranche ${sep=' -tranche ' recalibration_tranche_values} \
      -an ${sep=' -an ' recalibration_annotation_values} \
      -mode SNP \
      --input-model ${model_report} \
      -output-tranches-for-scatter \
      --max-gaussians 6 \
      -resource hapmap,known=false,training=true,truth=true,prior=15:${hapmap_resource_vcf} \
      -resource omni,known=false,training=true,truth=true,prior=12:${omni_resource_vcf} \
      -resource 1000G,known=false,training=true,truth=false,prior=10:${one_thousand_genomes_resource_vcf} \
      -resource dbsnp,known=true,training=false,truth=false,prior=7:${dbsnp_resource_vcf}
  }
  runtime {
    docker: docker_image
    memory: mem_size
    cpu: "2"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: preemptibles
  }
  output {
    File recalibration = "${recalibration_filename}"
    File recalibration_index = "${recalibration_filename}.idx"
    File tranches = "${tranches_filename}"
  }
}

task GatherTranches {
  File input_fofn
  String output_filename

  String gatk_path
  String java_opt

  String docker_image
  Int disk_size
  String mem_size
  Int preemptibles

  command <<<
    set -e
    set -o pipefail

    # this is here to deal with the JES bug where commands may be run twice
    rm -rf tranches

    mkdir tranches
    RETRY_LIMIT=5

    count=0
    until cat ${input_fofn} | /google-cloud-sdk/bin/gsutil -m cp -L cp.log -c -I tranches/; do
        sleep 1
        ((count++)) && ((count >= $RETRY_LIMIT)) && break
    done
    if [ "$count" -ge "$RETRY_LIMIT" ]; then
        echo 'Could not copy all the tranches from the cloud' && exit 1
    fi

    #cat ${input_fofn} | rev | cut -d '/' -f 1 | rev | awk '{print "tranches/" $1}' > inputs.list #.list has changed to .args for gatk4, but will change back to .list soon
    cat ${input_fofn} | rev | cut -d '/' -f 1 | rev | awk '{print "tranches/" $1}' > inputs.args

      ${gatk_path} --java-options "${java_opt}" \
      GatherTranches \
      --input inputs.args \
      --output ${output_filename}
  >>>
  runtime {
    docker: docker_image
    memory: mem_size
    cpu: "2"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: preemptibles
  }
  output {
    File tranches = "${output_filename}"
  }
}

task ApplyRecalibration {
  String recalibrated_vcf_filename
  File input_vcf
  File input_vcf_index
  File indels_recalibration
  File indels_recalibration_index
  File indels_tranches
  File snps_recalibration
  File snps_recalibration_index
  File snps_tranches

  Float indel_filter_level
  Float snp_filter_level

  String gatk_path
  String java_opt

  String docker_image
  Int disk_size
  String mem_size
  Int preemptibles

  command {
    set -e

    ${gatk_path} --java-options "${java_opt}" \
      ApplyVQSR \
      -O tmp.indel.recalibrated.vcf \
      -V ${input_vcf} \
      --recal-file ${indels_recalibration} \
      -tranches-file ${indels_tranches} \
      -truth-sensitivity-filter-level ${indel_filter_level} \
      --create-output-variant-index true \
      -mode INDEL

    ${gatk_path} --java-options "${java_opt}" \
      ApplyVQSR \
      -O ${recalibrated_vcf_filename} \
      -V tmp.indel.recalibrated.vcf \
      --recal-file ${snps_recalibration} \
      -tranches-file ${snps_tranches} \
      -truth-sensitivity-filter-level ${snp_filter_level} \
      --create-output-variant-index true \
      -mode SNP
  }
  runtime {
    docker: docker_image
    memory: mem_size
    cpu: "1"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: preemptibles
  }
  output {
    File recalibrated_vcf = "${recalibrated_vcf_filename}"
    File recalibrated_vcf_index = "${recalibrated_vcf_filename}.tbi"
  }
}

task GatherVcfs {
  File input_vcfs_fofn
  String output_vcf_name

  String gatk_path
  String java_opt

  String docker_image
  Int disk_size
  String mem_size
  Int preemptibles

  command <<<
    set -e

    # Now using NIO to localize the vcfs but the input file must have a ".list" extension
    #mv ${input_vcfs_fofn} inputs.list #.list has changed to .args for gatk4, but will change back to .list soon
    mv ${input_vcfs_fofn} inputs.args

    # ignoreSafetyChecks make a big performance difference so we include it in our invocation
    ${gatk_path} --java-options "${java_opt}" \
    GatherVcfsCloud \
    --ignore-safety-checks \
    --gather-type BLOCK \
    --input inputs.args \
    --output ${output_vcf_name}

    ${gatk_path} --java-options "-Xmx6g -Xms6g" \
    IndexFeatureFile \
    --feature-file ${output_vcf_name}
  >>>
  runtime {
    docker: docker_image
    memory: mem_size
    cpu: "1"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: preemptibles
  }
  output {
    File output_vcf = "${output_vcf_name}"
    File output_vcf_index = "${output_vcf_name}.tbi"
  }
}

task CollectVariantCallingMetrics {
  File input_vcf
  File input_vcf_index

  String metrics_filename_prefix
  File dbsnp_vcf
  File dbsnp_vcf_index
  File interval_list
  File ref_dict

  String gatk_path
  String java_opt

  String docker_image
  Int disk_size
  String mem_size
  Int preemptibles

  command {
    ${gatk_path} --java-options "${java_opt}" \
      CollectVariantCallingMetrics \
      --INPUT ${input_vcf} \
      --DBSNP ${dbsnp_vcf} \
      --SEQUENCE_DICTIONARY ${ref_dict} \
      --OUTPUT ${metrics_filename_prefix} \
      --THREAD_COUNT 8 \
      --TARGET_INTERVALS ${interval_list}
  }
  output {
    File detail_metrics_file = "${metrics_filename_prefix}.variant_calling_detail_metrics"
    File summary_metrics_file = "${metrics_filename_prefix}.variant_calling_summary_metrics"
  }
  runtime {
    docker: docker_image
    memory: mem_size
    cpu: 2
    disks: "local-disk " + disk_size + " HDD"
    preemptible: preemptibles
  }
}

task GatherMetrics {
  File input_details_fofn
  File input_summaries_fofn

  String output_prefix

  String gatk_path
  String java_opt

  String docker_image
  Int disk_size
  String mem_size
  Int preemptibles

  command <<<
    set -e
    set -o pipefail

    # this is here to deal with the JES bug where commands may be run twice
    rm -rf metrics

    mkdir metrics
    RETRY_LIMIT=5

    count=0
    until cat ${input_details_fofn} | /google-cloud-sdk/bin/gsutil -m cp -L cp.log -c -I metrics/; do
        sleep 1
        ((count++)) && ((count >= $RETRY_LIMIT)) && break
    done
    if [ "$count" -ge "$RETRY_LIMIT" ]; then
        echo 'Could not copy all the metrics from the cloud' && exit 1
    fi

    count=0
    until cat ${input_summaries_fofn} | /google-cloud-sdk/bin/gsutil -m cp -L cp.log -c -I metrics/; do
        sleep 1
        ((count++)) && ((count >= $RETRY_LIMIT)) && break
    done
    if [ "$count" -ge "$RETRY_LIMIT" ]; then
        echo 'Could not copy all the metrics from the cloud' && exit 1
    fi

    INPUT=`cat ${input_details_fofn} | rev | cut -d '/' -f 1 | rev | sed s/.variant_calling_detail_metrics//g | awk '{printf("I=metrics/%s ", $1)}'`

    ${gatk_path} --java-option "${java_opt}" \
    AccumulateVariantCallingMetrics \
    --INPUT $INPUT \
    --OUTPUT ${output_prefix}
  >>>
  runtime {
    docker: docker_image
    memory: mem_size
    cpu: "1"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: preemptibles
  }
  output {
    File detail_metrics_file = "${output_prefix}.variant_calling_detail_metrics"
    File summary_metrics_file = "${output_prefix}.variant_calling_summary_metrics"
  }
}

task DynamicallyCombineIntervals {
  File intervals
  Int merge_count
  String docker_image
  String mem_size
  Int preemptibles

  command {
    python << CODE
    def parse_interval(interval):
        colon_split = interval.split(":")
        chromosome = colon_split[0]
        dash_split = colon_split[1].split("-")
        start = int(dash_split[0])
        end = int(dash_split[1])
        return chromosome, start, end

    def add_interval(chr, start, end):
        lines_to_write.append(chr + ":" + str(start) + "-" + str(end))
        return chr, start, end

    count = 0
    chain_count = ${merge_count}
    l_chr, l_start, l_end = "", 0, 0
    lines_to_write = []
    with open("${intervals}") as f:
        with open("out.intervals", "w") as f1:
            for line in f.readlines():
                # initialization
                if count == 0:
                    w_chr, w_start, w_end = parse_interval(line)
                    count = 1
                    continue
                # reached number to combine, so spit out and start over
                if count == chain_count:
                    l_char, l_start, l_end = add_interval(w_chr, w_start, w_end)
                    w_chr, w_start, w_end = parse_interval(line)
                    count = 1
                    continue

                c_chr, c_start, c_end = parse_interval(line)
                # if adjacent keep the chain going
                if c_chr == w_chr and c_start == w_end + 1:
                    w_end = c_end
                    count += 1
                    continue
                # not adjacent, end here and start a new chain
                else:
                    l_char, l_start, l_end = add_interval(w_chr, w_start, w_end)
                    w_chr, w_start, w_end = parse_interval(line)
                    count = 1
            if l_char != w_chr or l_start != w_start or l_end != w_end:
                add_interval(w_chr, w_start, w_end)
            f1.writelines("\n".join(lines_to_write))
    CODE
  }

  runtime {
    memory: mem_size
    preemptible: preemptibles
    docker: docker_image
  }

  output {
    File output_intervals = "out.intervals"
  }
}
