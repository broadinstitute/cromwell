
workflow Arrays {
  String sample_alias
  Int analysis_version_number
  Float call_rate_threshold
  String reported_gender

  String idat_dir_name
  File file_of_idat_filenames
  File ref_fasta
  File ref_fasta_index
  File ref_dict

  File dbSNP_vcf
  File dbSNP_vcf_index

  File params_file

  File bead_pool_manifest_file
  String bead_pool_manifest_filename = sub(bead_pool_manifest_file, "gs://.*/", "")
  String chip_type = sub(bead_pool_manifest_filename, "\\.bpm$", "")

  File extended_chip_manifest_file
  File cluster_file
  File? gender_cluster_file
  File? zcall_thresholds_file

  # For CheckFingerprint:
  File? fingerprint_genotypes_vcf_file
  File? fingerprint_genotypes_vcf_index_file
  File haplotype_database_file

  # For SelectVariants
  File variant_rsids_file

  # For HapMap GenotypeConcordance Check:
  File? control_sample_vcf_file
  File? control_sample_vcf_index_file
  File? control_sample_intervals_file
  String? control_sample_name

  Int disk_size
  Int preemptible_tries

  call UpdateChipWellBarcodeIndex {
    input:
      params_file = params_file,
      disk_size = disk_size,
      preemptible_tries = preemptible_tries
  }

  call AutoCall {
    input:
      idat_dir_name = idat_dir_name,
      file_of_idat_filenames = file_of_idat_filenames,
      bead_pool_manifest_file = bead_pool_manifest_file,
      bead_pool_manifest_filename = bead_pool_manifest_filename,
      cluster_file = cluster_file,
      disk_size = disk_size,
      preemptible_tries = preemptible_tries
  }

  if (size(AutoCall.gtc_file) == 0) {
    call GetFailedAutocallVersion {
      input:
        autocallStdOut = AutoCall.autocall_stdout
    }

    call GenerateEmptyVariantCallingMetricsFile {
      input:
        chip_well_barcode = idat_dir_name,
        sample_alias = sample_alias,
        chip_type = chip_type,
        reported_gender = reported_gender,
        autocall_version = GetFailedAutocallVersion.autocallVersion,
        output_metrics_basename = sample_alias,
        cluster_file = cluster_file,
        analysis_version_number = analysis_version_number,
        preemptible_tries = preemptible_tries
    }

    call UploadArraysMetrics as UploadEmptyArraysMetrics {
      input:
        arrays_variant_calling_detail_metrics = GenerateEmptyVariantCallingMetricsFile.detail_metrics,
        disk_size = disk_size,
        preemptible_tries = preemptible_tries
    }

    call BlacklistBarcode {
      input:
        upload_metrics_output = UploadEmptyArraysMetrics.upload_metrics_empty_file,
        chip_well_barcode = idat_dir_name,
        preemptible_tries = preemptible_tries
    }
  }

  if (defined(gender_cluster_file)) {
      call AutoCall as GenderAutocall {
        input:
          idat_dir_name = idat_dir_name,
          file_of_idat_filenames = file_of_idat_filenames,
          bead_pool_manifest_file = bead_pool_manifest_file,
          bead_pool_manifest_filename = bead_pool_manifest_filename,
          cluster_file = gender_cluster_file,
          disk_size = disk_size,
          preemptible_tries = preemptible_tries
      }
  }

 if (size(AutoCall.gtc_file) > 0) {
  call GtcToVcf {
    input:
      vcf_filename = sub(sub(AutoCall.gtc_file, "gs://.*/", ""), ".gtc$", "") + ".vcf.gz",
      input_gtc = AutoCall.gtc_file,
      gender_gtc = GenderAutocall.gtc_file,
      extended_chip_manifest_file = extended_chip_manifest_file,
      cluster_file = cluster_file,
      normalization_manifest_file = AutoCall.bead_pool_manifest_csv_file,
      zcall_thresholds_file = zcall_thresholds_file,
      sample_alias = sample_alias,
      analysis_version_number = analysis_version_number,
      reported_gender = reported_gender,
      fingerprint_genotypes_vcf_file = fingerprint_genotypes_vcf_file,
      fingerprint_genotypes_vcf_index_file = fingerprint_genotypes_vcf_index_file,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      disk_size = disk_size,
      preemptible_tries = preemptible_tries
  }

  if (defined(zcall_thresholds_file)) {
    call zCall {
      input:
        zcall_ped_filename = sub(sub(AutoCall.gtc_file, "gs://.*/", ""), ".gtc$", "") + ".zcall.ped",
        zcall_map_filename = sub(sub(AutoCall.gtc_file, "gs://.*/", ""), ".gtc$", "") + ".zcall.map",
        input_gtc = AutoCall.gtc_file,
        bead_pool_manifest_csv_file = AutoCall.bead_pool_manifest_csv_file,
        zcall_thresholds_file = zcall_thresholds_file,
        disk_size = disk_size,
        preemptible_tries = preemptible_tries
    }
  }

  if (defined(zCall.ped_file) && defined(zCall.map_file)) {
    call MergePedIntoVcf {
      input:
        input_vcf = GtcToVcf.output_vcf,
        input_vcf_index = GtcToVcf.output_vcf_index,
        output_vcf_filename = sub(GtcToVcf.output_vcf, "gs://.*/", ""),
        ped_file = zCall.ped_file,
        map_file = zCall.map_file,
        disk_size = disk_size,
        preemptible_tries = preemptible_tries
    }
  }

  # if zCall doesn't run, then MergePedIntoVcf doesn't run.
  # if MergePedIntoVcf doesn't run, then MergePedIntoVcf.output_vcf and output_vcf_index should just be
  # the MergePedIntoVcf.input_vcf and input_vcf_index
  # using select_first to cast File? to File type
  File MergePedIntoVcfOutputVcf = if (defined(zCall.ped_file) && defined(zCall.map_file)) then select_first([MergePedIntoVcf.output_vcf]) else GtcToVcf.output_vcf
  File MergePedIntoVcfOutputVcfIndex = if (defined(zCall.ped_file) && defined(zCall.map_file)) then select_first([MergePedIntoVcf.output_vcf_index]) else GtcToVcf.output_vcf_index

  call CollectArraysVariantCallingMetrics {
    input:
      input_vcf_file = MergePedIntoVcfOutputVcf,
      input_vcf_index_file = MergePedIntoVcfOutputVcfIndex,
      dbSNP_vcf_file = dbSNP_vcf,
      dbSNP_vcf_index_file = dbSNP_vcf_index,
      call_rate_threshold = call_rate_threshold,
      output_metrics_basename = sub(sub(AutoCall.gtc_file, "gs://.*/", ""), ".gtc$", ""),
      disk_size = disk_size,
      preemptible_tries = preemptible_tries
  }

  call VcfToIntervalList {
    input:
      vcf_file = MergePedIntoVcfOutputVcf,
      interval_list_file_name = sub(sub(MergePedIntoVcfOutputVcf, "gs://.*/", ""), ".vcf.gz$", "") + ".interval_list",
      disk_size = disk_size,
      preemptible_tries = preemptible_tries
  }

  if (defined(fingerprint_genotypes_vcf_file) && defined(fingerprint_genotypes_vcf_index_file)) {
    call CheckFingerprint {
      input:
        input_vcf_file = MergePedIntoVcfOutputVcf,
        input_vcf_index_file = MergePedIntoVcfOutputVcfIndex,
        genotypes_vcf_file = fingerprint_genotypes_vcf_file,
        genotypes_vcf_index_file = fingerprint_genotypes_vcf_index_file,
        haplotype_database_file = haplotype_database_file,
        observed_sample_alias = sub(sub(AutoCall.gtc_file, "gs://.*/", ""), ".gtc$", ""),
        expected_sample_alias = sample_alias,
        output_metrics_basename = sub(sub(AutoCall.gtc_file, "gs://.*/", ""), ".gtc$", ""),
        disk_size = disk_size,
        preemptible_tries = preemptible_tries
    }

    call SelectVariants {
      input:
        input_vcf_file = MergePedIntoVcfOutputVcf,
        input_vcf_index_file = MergePedIntoVcfOutputVcfIndex,
        variant_rsids_file = variant_rsids_file,
        disk_size = disk_size,
        preemptible_tries = preemptible_tries
    }
  }

  if (defined(control_sample_vcf_file) && defined(control_sample_vcf_index_file) &&
      defined(control_sample_intervals_file) && defined(control_sample_name)) {
    call GenotypeConcordance {
      input:
        call_vcf_file = MergePedIntoVcfOutputVcf,
        call_vcf_index_file = MergePedIntoVcfOutputVcfIndex,
        call_intervals_file = VcfToIntervalList.interval_list_file,
        call_sample_name = sub(sub(AutoCall.gtc_file, "gs://.*/", ""), ".gtc$", ""),
        truth_vcf_file = control_sample_vcf_file,
        truth_vcf_index_file = control_sample_vcf_index_file,
        truth_intervals_file = control_sample_intervals_file,
        truth_sample_name = control_sample_name,
        output_metrics_basename = sub(sub(AutoCall.gtc_file, "gs://.*/", ""), ".gtc$", ""),
        disk_size = disk_size,
        preemptible_tries = preemptible_tries
    }
  }

  call UploadArraysMetrics {
    input:
      arrays_variant_calling_detail_metrics = CollectArraysVariantCallingMetrics.detail_metrics,
      arrays_variant_calling_summary_metrics = CollectArraysVariantCallingMetrics.summary_metrics,
      arrays_control_code_summary_metrics = CollectArraysVariantCallingMetrics.control_metrics,
      fingerprinting_detail_metrics = CheckFingerprint.detail_metrics,
      fingerprinting_summary_metrics = CheckFingerprint.summary_metrics,
      genotype_concordance_summary_metrics = GenotypeConcordance.summary_metrics,
      genotype_concordance_detail_metrics  = GenotypeConcordance.detail_metrics,
      genotype_concordance_contingency_metrics = GenotypeConcordance.contingency_metrics,
      disk_size = disk_size,
      preemptible_tries = preemptible_tries
  }
}
  File ArraysVariantCallingDetailMetrics = if (defined(CollectArraysVariantCallingMetrics.detail_metrics)) then select_first([CollectArraysVariantCallingMetrics.detail_metrics]) else select_first([GenerateEmptyVariantCallingMetricsFile.detail_metrics])

  output {
    AutoCall.gtc_file
    File ArraysVariantCallingDetailMetricsFile = ArraysVariantCallingDetailMetrics
    File? ArraysVariantCallingSummaryMetricsFile = CollectArraysVariantCallingMetrics.summary_metrics
    File? ArraysVariantCallingControlMetricsFile = CollectArraysVariantCallingMetrics.control_metrics
    File?? FingerprintDetailMetricsFile = CheckFingerprint.detail_metrics
    File?? FingerprintSummaryMetricsFile = CheckFingerprint.summary_metrics
    File?? GenotypeConcordanceSummaryMetricsFile = GenotypeConcordance.summary_metrics
    File?? GenotypeConcordanceDetailMetricsFile  = GenotypeConcordance.detail_metrics
    File?? GenotypeConcordanceContingencyMetricsFile = GenotypeConcordance.contingency_metrics
    File? MergePedIntoVcfOutputVcfFile = MergePedIntoVcfOutputVcf
    File? MergePedIntoVcfOutputVcfIndexFile = MergePedIntoVcfOutputVcfIndex
  }
}

task AutoCall {
  String idat_dir_name
  File file_of_idat_filenames
  File bead_pool_manifest_file
  String bead_pool_manifest_filename
  File? cluster_file
  Int disk_size
  Int preemptible_tries

  command <<<
    set -e
    rm -rf ${idat_dir_name}
    mkdir ${idat_dir_name}

    RETRY_LIMIT=5

    until cat ${file_of_idat_filenames} | /root/google-cloud-sdk/bin/gsutil -m cp -L cp.log -c -I ${idat_dir_name}/; do
        sleep 1
        ((count++)) && ((count==$RETRY_LIMIT)) && break
    done

    if [ "$count" = "$RETRY_LIMIT" ]; then
        echo 'Could not copy all the files from the cloud' && exit 1
    fi

    mono /usr/gitc/autocall/AutoConvert.exe ${idat_dir_name} ${idat_dir_name} ${bead_pool_manifest_file} ${cluster_file} | tee autocall_std.txt

    #wtf autocall writes the manifest csv to the root directory with an escape character before it e.g. /\PsychChip_v1-1_15073391_A1.bpm.csv  WHY!!
    cp /*.bpm.csv ${bead_pool_manifest_filename}.csv

    if grep -q "Normalization failed" autocall_std.txt; then
      # make an empty gtc file so that jes won't fail this task and wf can move on
      touch ${idat_dir_name}/${idat_dir_name}.gtc
    fi
  >>>
  runtime {
    disks: "local-disk " + disk_size + " HDD"
    memory: "10 GB"
    preemptible: preemptible_tries
  }
  output {
    File gtc_file = "${idat_dir_name}/${idat_dir_name}.gtc"
    File bead_pool_manifest_csv_file = "${bead_pool_manifest_filename}.csv"
    File autocall_stdout = "autocall_std.txt"
  }
}

task GetFailedAutocallVersion {
  File autocallStdOut

  command <<<
    stringContainsVersion=`grep "Calling Utility" ${autocallStdOut}`
    echo "$stringContainsVersion" | cut -d ' ' -f 5-
  >>>

  output {
    String autocallVersion = read_string(stdout())
  }
}

task GenerateEmptyVariantCallingMetricsFile {
  String chip_well_barcode
  String sample_alias
  String chip_type
  String reported_gender
  String autocall_version
  String output_metrics_basename
  File? cluster_file
  Int analysis_version_number
  Int preemptible_tries

  command <<<
    java -Xmx2g -jar /usr/gitc/picard-private.jar \
            GenerateEmptyVariantCallingMetrics \
            CHIP_WELL_BARCODE=${chip_well_barcode} \
            SAMPLE_ALIAS="${sample_alias}" \
            CHIP_TYPE=${chip_type} \
            REPORTED_GENDER="${reported_gender}" \
            CLUSTER_FILE_NAME=${cluster_file} \
            AUTOCALL_VERSION=${autocall_version} \
            ANALYSIS_VERSION=${analysis_version_number} \
            OUTPUT=${output_metrics_basename}
  >>>

  runtime {
      memory: "2 GB"
      preemptible: preemptible_tries
  }

  output {
      File detail_metrics = "${output_metrics_basename}.arrays_variant_calling_detail_metrics"
    }
}

task BlacklistBarcode {
  File upload_metrics_output
  String chip_well_barcode
  Int preemptible_tries

  command <<<
    java -Xmx2g -jar /usr/gitc/picard-private.jar \
                  ArraysManualBlacklistUpdate \
                  CHIP_WELL_BARCODE=${chip_well_barcode} \
                  REASON=DATA_QUALITY \
                  DB_USERNAME_FILE=/usr/gitc/cloudsql.db_username.txt \
                  DB_PASSWORD_FILE=/usr/gitc/cloudsql.db_password.txt \
                  DB_JDBC_FILE=/usr/gitc/cloudsql.db_jdbc.txt \
                  NOTES="Normalization failed"
  >>>
  runtime {
    memory: "2 GB"
    preemptible: preemptible_tries
  }
}

task zCall {
  String zcall_ped_filename
  String zcall_map_filename
  File input_gtc
  File bead_pool_manifest_csv_file
  File? zcall_thresholds_file
  Int disk_size
  Int preemptible_tries


  command <<<
    python /usr/gitc/zcall/zCall.py -B ${bead_pool_manifest_csv_file} -G ${input_gtc} -T ${zcall_thresholds_file} > ${zcall_ped_filename}
    python /usr/gitc/zcall/makeMAPfile.py -B ${bead_pool_manifest_csv_file} > ${zcall_map_filename}
  >>>
  runtime {
    disks: "local-disk " + disk_size + " HDD"
    memory: "10 GB"
    preemptible: preemptible_tries
  }
  output {
    File? ped_file = "${zcall_ped_filename}"
    File? map_file = "${zcall_map_filename}"
  }
}

task GtcToVcf {
  String vcf_filename
  File input_gtc
  File? gender_gtc
  File extended_chip_manifest_file
  File cluster_file
  File normalization_manifest_file
  File? zcall_thresholds_file
  String sample_alias
  Int analysis_version_number
  String reported_gender

  File ref_fasta
  File ref_fasta_index
  File ref_dict

  Int disk_size
  Int preemptible_tries

  File? fingerprint_genotypes_vcf_file
  File? fingerprint_genotypes_vcf_index_file

  command {
  java -Xmx7g -jar /usr/gitc/picard-private.jar \
    GtcToVcf \
    INPUT=${input_gtc} \
    ${"GENDER_GTC=" + gender_gtc} \
    ${"ZCALL_THRESHOLDS_FILE=" + zcall_thresholds_file} \
    OUTPUT=${vcf_filename} \
    MANIFEST=${extended_chip_manifest_file} \
    CLUSTER_FILE=${cluster_file} \
    ILLUMINA_NORMALIZATION_MANIFEST=${normalization_manifest_file} \
    SAMPLE_ALIAS="${sample_alias}" \
    ANALYSIS_VERSION_NUMBER=${analysis_version_number} \
    EXPECTED_GENDER="${reported_gender}" \
    ${"FINGERPRINT_GENOTYPES_VCF_FILE=" + fingerprint_genotypes_vcf_file} \
    REFERENCE_SEQUENCE=${ref_fasta} \
    CREATE_INDEX=true
  }
  runtime {
    disks: "local-disk " + disk_size + " HDD"
    memory: "7 GB"
    preemptible: preemptible_tries
  }
  output {
    File output_vcf = "${vcf_filename}"
    File output_vcf_index = "${vcf_filename}.tbi"
  }
}

task CollectArraysVariantCallingMetrics {
  File input_vcf_file
  File input_vcf_index_file
  File dbSNP_vcf_file
  File dbSNP_vcf_index_file
  Float call_rate_threshold
  String output_metrics_basename

  Int disk_size
  Int preemptible_tries

  command {
    java -Xmx2g -jar /usr/gitc/picard-private.jar \
      CollectArraysVariantCallingMetrics \
      INPUT=${input_vcf_file} \
      DBSNP=${dbSNP_vcf_file} \
      CALL_RATE_PF_THRESHOLD=${call_rate_threshold} \
      OUTPUT=${output_metrics_basename}
  }
  runtime {
    disks: "local-disk " + disk_size + " HDD"
    memory: "2 GB"
    preemptible: preemptible_tries
  }
  output {
    File summary_metrics = "${output_metrics_basename}.arrays_variant_calling_summary_metrics"
    File detail_metrics = "${output_metrics_basename}.arrays_variant_calling_detail_metrics"
    File control_metrics = "${output_metrics_basename}.arrays_control_code_summary_metrics"
  }
}

task VcfToIntervalList {
  File vcf_file
  String interval_list_file_name

  Int disk_size
  Int preemptible_tries

  command {
    java -Xmx2g -jar /usr/gitc/picard-private.jar \
      VcfToIntervalList \
      INPUT=${vcf_file} \
      OUTPUT=${interval_list_file_name}
  }
  runtime {
    disks: "local-disk " + disk_size + " HDD"
    memory: "2 GB"
    preemptible: preemptible_tries
  }
  output {
    File interval_list_file = "${interval_list_file_name}"
  }
}

task CheckFingerprint {
  File input_vcf_file
  File input_vcf_index_file
  File? genotypes_vcf_file
  File? genotypes_vcf_index_file
  File haplotype_database_file
  String observed_sample_alias
  String expected_sample_alias
  String output_metrics_basename

  Int disk_size
  Int preemptible_tries

  # Paraphrased from Yossi:
  # Override the default LOD threshold of 5 because if the PL field
  # is missing from the VCF, CheckFingerprint will default to an error
  # rate equivalent to a LOD score of 2, and we don't want to see
  # confident LOD scores w/ no confident SNPs.
  Float genotype_lod_threshold = 1.9

  command <<<
    java -Xmx2g -jar /usr/gitc/picard-private.jar \
      CheckFingerprint \
      INPUT=${input_vcf_file} \
      OBSERVED_SAMPLE_ALIAS="${observed_sample_alias}" \
      ${"GENOTYPES=" + genotypes_vcf_file} \
      EXPECTED_SAMPLE_ALIAS="${expected_sample_alias}" \
      HAPLOTYPE_MAP=${haplotype_database_file} \
      GENOTYPE_LOD_THRESHOLD=${genotype_lod_threshold} \
      OUTPUT=${output_metrics_basename}
  >>>
 runtime {
    disks: "local-disk " + disk_size + " HDD"
    memory: "2 GB"
    preemptible: preemptible_tries
  }
  output {
    File summary_metrics = "${output_metrics_basename}.fingerprinting_summary_metrics"
    File detail_metrics = "${output_metrics_basename}.fingerprinting_detail_metrics"
  }
}

task SelectVariants {
  File input_vcf_file
  File input_vcf_index_file
  File variant_rsids_file

  Int disk_size
  Int preemptible_tries

  String fingerprints_file = "fingerprints.vcf.gz"

  command <<<
    export GATK_LOCAL_JAR="/root/gatk.jar"
    gatk --java-options "-Xms2g" \
      SelectVariants \
      -V ${input_vcf_file} \
      --keep-ids ${variant_rsids_file} \
      -O ${fingerprints_file}
  >>>
  runtime {
    disks: "local-disk " + disk_size + " HDD"
    docker: "us.gcr.io/broad-gatk/gatk:4.0.0.0"
    memory: "3.5 GB"
    preemptible: preemptible_tries
  }
  output {
    File fingerprints = fingerprints_file
    File fingerprints_index = fingerprints_file + ".tbi"
  }
}

task GenotypeConcordance {
  File call_vcf_file
  File call_vcf_index_file
  File call_intervals_file
  String call_sample_name
  File? truth_vcf_file
  File? truth_vcf_index_file
  File? truth_intervals_file
  String? truth_sample_name
  String output_metrics_basename

  Int disk_size
  Int preemptible_tries

  command {
    java -Xmx2g -jar /usr/gitc/picard-private.jar \
      GenotypeConcordance \
      CALL_VCF=${call_vcf_file} \
      CALL_SAMPLE=${call_sample_name} \
      TRUTH_VCF=${truth_vcf_file} \
      TRUTH_SAMPLE=${truth_sample_name} \
      INTERVALS=${call_intervals_file} \
      INTERVALS=${truth_intervals_file} \
      OUTPUT=${output_metrics_basename}
  }
  runtime {
    disks: "local-disk " + disk_size + " HDD"
    memory: "2 GB"
    preemptible: preemptible_tries
  }
  output {
    File summary_metrics = "${output_metrics_basename}.genotype_concordance_summary_metrics"
    File detail_metrics = "${output_metrics_basename}.genotype_concordance_detail_metrics"
    File contingency_metrics = "${output_metrics_basename}.genotype_concordance_contingency_metrics"
  }
}

task UploadArraysMetrics {
  File arrays_variant_calling_detail_metrics
  File? arrays_variant_calling_summary_metrics
  File? arrays_control_code_summary_metrics
  File? fingerprinting_detail_metrics
  File? fingerprinting_summary_metrics
  File? genotype_concordance_summary_metrics
  File? genotype_concordance_detail_metrics
  File? genotype_concordance_contingency_metrics

  Int disk_size
  Int preemptible_tries

  command <<<
    rm -rf metrics_upload_dir &&
    mkdir metrics_upload_dir &&

    # check that files are passed in before copying them -- [ -z FILE ] evaluates to true if FILE not there
    ! [ -z ${genotype_concordance_summary_metrics} ] &&
    cp ${genotype_concordance_summary_metrics} metrics_upload_dir
    ! [ -z ${genotype_concordance_detail_metrics} ] &&
    cp ${genotype_concordance_detail_metrics} metrics_upload_dir
    ! [ -z ${genotype_concordance_contingency_metrics} ] &&
    cp ${genotype_concordance_contingency_metrics} metrics_upload_dir

    ! [ -z ${fingerprinting_detail_metrics} ] &&
    cp ${fingerprinting_detail_metrics} metrics_upload_dir
    ! [ -z ${fingerprinting_summary_metrics} ] &&
    cp ${fingerprinting_summary_metrics} metrics_upload_dir

    cp ${arrays_variant_calling_detail_metrics} metrics_upload_dir
    ! [ -z ${arrays_variant_calling_summary_metrics} ] &&
    cp ${arrays_variant_calling_summary_metrics} metrics_upload_dir

    ! [ -z ${arrays_control_code_summary_metrics} ] &&
    cp ${arrays_control_code_summary_metrics} metrics_upload_dir
    java -Xmx2g -jar /usr/gitc/picard-private.jar \
      UploadArraysMetrics \
      ANALYSIS_DIRECTORY=metrics_upload_dir \
      DB_USERNAME_FILE=/usr/gitc/cloudsql.db_username.txt \
      DB_PASSWORD_FILE=/usr/gitc/cloudsql.db_password.txt \
      DB_JDBC_FILE=/usr/gitc/cloudsql.db_jdbc.txt &&
    touch empty_file_for_dependency
  >>>
  runtime {
    disks: "local-disk " + disk_size + " HDD"
    memory: "2 GB"
    preemptible: preemptible_tries
  }
  output {
      File upload_metrics_empty_file = "empty_file_for_dependency"
    }
}

task UpdateChipWellBarcodeIndex {
  File params_file
  Int disk_size
  Int preemptible_tries

  command <<<
    java -Xmx2g -jar /usr/gitc/picard-private.jar \
      UpdateChipWellBarcodeIndex \
      PARAMS_FILE=${params_file} \
      DB_USERNAME_FILE=/usr/gitc/cloudsql.db_username.txt \
      DB_PASSWORD_FILE=/usr/gitc/cloudsql.db_password.txt \
      DB_JDBC_FILE=/usr/gitc/cloudsql.db_jdbc.txt

  >>>
  runtime {
    disks: "local-disk " + disk_size + " HDD"
    memory: "2 GB"
    preemptible: preemptible_tries
  }
}

task MergePedIntoVcf {
  File input_vcf
  File input_vcf_index
  File? ped_file
  File? map_file

  String output_vcf_filename

  Int disk_size
  Int preemptible_tries

  command {
    java -Xmx3g -jar /usr/gitc/picard-private.jar \
      MergePedIntoVcf \
      ORIGINAL_VCF=${input_vcf} \
      PED_FILE=${ped_file} \
      MAP_FILE=${map_file} \
      OUTPUT=${output_vcf_filename} \
      CREATE_INDEX=true
  }

  runtime {
    memory: "3500 MB"
    cpu: "1"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: preemptible_tries
  }
  output {
    File output_vcf = "${output_vcf_filename}"
    File output_vcf_index = "${output_vcf_filename}.tbi"
  }
}