task AutoCall {
  String idat_dir_name
  File file_of_idat_filenames
  File bead_pool_manifest_file
  String bead_pool_manifest_filename
  File cluster_file
  Int disk_size
  Int preemptible_tries

  command <<<

    rm -rf ${idat_dir_name}
    mkdir ${idat_dir_name}

    if [ -s ${cluster_file} ]; then
      RETRY_LIMIT=5

      until cat ${file_of_idat_filenames} | /root/google-cloud-sdk/bin/gsutil -m cp -L cp.log -c -I ${idat_dir_name}/; do
          sleep 1
          ((count++)) && ((count==$RETRY_LIMIT)) && break
      done

      if [ "$count" = "$RETRY_LIMIT" ]; then
          echo 'Could not copy all the files from the cloud' && exit 1
      fi

      mono /usr/gitc/autocall/AutoConvert.exe ${idat_dir_name} ${idat_dir_name} ${bead_pool_manifest_file} ${cluster_file}
      #wtf autocall writes the manifest csv to the root directory with an escape character before it e.g. /\PsychChip_v1-1_15073391_A1.bpm.csv  WHY!!
      cp /*.bpm.csv ${bead_pool_manifest_filename}.csv
    else
      echo "Found 0-length cluster_file, not running autocall"
      touch ${idat_dir_name}/${idat_dir_name}.gtc
      touch ${bead_pool_manifest_filename}.csv
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
  }
}

task zCall {
  String zcall_ped_filename
  String zcall_map_filename
  File input_gtc
  File bead_pool_manifest_csv_file
  File zcall_thresholds_file   # if this file is empty (0-length) the workflow will not run zCall
  Int disk_size
  Int preemptible_tries


  command <<<
    if [ -s ${zcall_thresholds_file} ]; then
      python /usr/gitc/zcall/zCall.py -B ${bead_pool_manifest_csv_file} -G ${input_gtc} -T ${zcall_thresholds_file} > ${zcall_ped_filename}
      python /usr/gitc/zcall/makeMAPfile.py -B ${bead_pool_manifest_csv_file} > ${zcall_map_filename}
    else
      echo "Found 0-length zCall thresholds file.  Not running zCall"
      touch ${zcall_ped_filename}
      touch ${zcall_map_filename}
    fi
  >>>
  runtime {
    disks: "local-disk " + disk_size + " HDD"
    memory: "10 GB"
    preemptible: preemptible_tries
  }
  output {
    File ped_file = "${zcall_ped_filename}"
    File map_file = "${zcall_map_filename}"
  }
}

task GtcToVcf {
  String vcf_filename
  File input_gtc
  File gender_gtc
  File extended_chip_manifest_file
  File cluster_file
  File normalization_manifest_file
  File zcall_thresholds_file
  String sample_alias
  Int analysis_version_number
  String reported_gender
  File fingerprint_gender_file

  File ref_fasta
  File ref_fasta_index
  File ref_dict

  Int disk_size
  Int preemptible_tries

  command {
  if [ -s ${gender_gtc} ]; then
    gender_gtc_file='GENDER_GTC=${gender_gtc}'
  fi

  if [ -s ${zcall_thresholds_file} ]; then
    zcall_file='ZCALL_THRESHOLDS_FILE=${zcall_thresholds_file}'
  fi

  java -Xmx7g -jar /usr/gitc/picard-private.jar \
    GtcToVcf \
    INPUT=${input_gtc} \
    $gender_gtc_file \
    $zcall_file \
    OUTPUT=${vcf_filename} \
    MANIFEST=${extended_chip_manifest_file} \
    CLUSTER_FILE=${cluster_file} \
    ILLUMINA_NORMALIZATION_MANIFEST=${normalization_manifest_file} \
    SAMPLE_ALIAS=${sample_alias} \
    ANALYSIS_VERSION_NUMBER=${analysis_version_number} \
    EXPECTED_GENDER=${reported_gender} \
    FINGERPRINT_GENDER_FILE=${fingerprint_gender_file} \
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
  File genotypes_vcf_file
  File genotypes_vcf_index_file
  File haplotype_database_file
  String observed_sample_alias
  String expected_sample_alias
  String output_metrics_basename

  Int disk_size
  Int preemptible_tries

  command <<<
  if [ -s ${genotypes_vcf_file} ]; then
    java -Xmx2g -jar /usr/gitc/picard-private.jar \
    CheckFingerprint \
    INPUT=${input_vcf_file} \
    OBSERVED_SAMPLE_ALIAS=${observed_sample_alias} \
    GENOTYPES=${genotypes_vcf_file} \
    EXPECTED_SAMPLE_ALIAS=${expected_sample_alias} \
    HAPLOTYPE_MAP=${haplotype_database_file} \
    OUTPUT=${output_metrics_basename}
  else
    echo "No fingerprint found. Skipping Fingerprint check"
    # We touch the outputs here in order to create 0 length files.  Otherwise the task will fail since the expected outputs are not to be found
    touch ${output_metrics_basename}.fingerprinting_summary_metrics
    touch ${output_metrics_basename}.fingerprinting_detail_metrics
  fi
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

task GenotypeConcordance {
  File call_vcf_file
  File call_vcf_index_file
  File call_intervals_file
  String call_sample_name
  File truth_vcf_file
  File truth_vcf_index_file
  File truth_intervals_file
  String truth_sample_name
  String output_metrics_basename

  Int disk_size
  Int preemptible_tries

  command {
    if [ -s ${truth_vcf_file} ]; then

      java -Xmx2g -jar /usr/gitc/picard-private.jar \
        GenotypeConcordance \
        CALL_VCF=${call_vcf_file} \
        CALL_SAMPLE=${call_sample_name} \
        TRUTH_VCF=${truth_vcf_file} \
        TRUTH_SAMPLE=${truth_sample_name} \
        INTERVALS=${call_intervals_file} \
        INTERVALS=${truth_intervals_file} \
        OUTPUT=${output_metrics_basename}

    else
      echo "No truth_sample_name provided. Skipping GenotypeConcordance"
      touch "${output_metrics_basename}.genotype_concordance_summary_metrics"
      touch "${output_metrics_basename}.genotype_concordance_detail_metrics"
      touch "${output_metrics_basename}.genotype_concordance_contingency_metrics"
    fi
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
  File arrays_variant_calling_summary_metrics
  File arrays_control_code_summary_metrics
  File fingerprinting_detail_metrics
  File fingerprinting_summary_metrics
  File genotype_concordance_summary_metrics
  File genotype_concordance_detail_metrics
  File genotype_concordance_contingency_metrics

  Int disk_size
  Int preemptible_tries

  command <<<
    rm -rf metrics_upload_dir &&
    mkdir metrics_upload_dir &&
    cp ${arrays_variant_calling_detail_metrics} metrics_upload_dir &&
    cp ${arrays_variant_calling_summary_metrics} metrics_upload_dir &&
    cp ${arrays_control_code_summary_metrics} metrics_upload_dir &&
    cp ${fingerprinting_detail_metrics} metrics_upload_dir &&
    cp ${fingerprinting_summary_metrics} metrics_upload_dir &&
    cp ${genotype_concordance_summary_metrics} metrics_upload_dir &&
    cp ${genotype_concordance_detail_metrics} metrics_upload_dir &&
    cp ${genotype_concordance_contingency_metrics} metrics_upload_dir &&
    java -Xmx2g -jar /usr/gitc/picard-private.jar \
      UploadArraysMetrics \
      ANALYSIS_DIRECTORY=metrics_upload_dir \
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
  File ped_file
  File map_file

  String output_vcf_filename

  Int disk_size
  Int preemptible_tries

  command {
    if [ -s ${ped_file} ]; then

      java -Xmx3g -jar /usr/gitc/picard-private.jar \
        MergePedIntoVcf \
        ORIGINAL_VCF=${input_vcf} \
        PED_FILE=${ped_file} \
        MAP_FILE=${map_file} \
        OUTPUT=${output_vcf_filename} \
        CREATE_INDEX=true

    else
      echo "0-length ped file found, not running MergePedIntoVcf"
      cp ${input_vcf} ${output_vcf_filename}
      cp ${input_vcf_index} ${output_vcf_filename}.tbi
    fi
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

  File extended_chip_manifest_file
  File cluster_file
  File gender_cluster_file
  File zcall_thresholds_file

  # For GtcToVcf
  File fingerprint_gender_file

  # For CheckFingerprint:
  File fingerprint_genotypes_vcf_file # if this file is empty (0-length) the workflow should not do fingerprint comparison (as there are no fingerprints for the sample)
  File fingerprint_genotypes_vcf_index_file
  File haplotype_database_file

  # For HapMap GenotypeConcordance Check:
  File control_sample_vcf_file
  File control_sample_vcf_index_file
  File control_sample_intervals_file
  String control_sample_name


  Int disk_size
  Int preemptible_tries

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

  call GtcToVcf {
    input:
      vcf_filename =  sub(sub(AutoCall.gtc_file, "gs://.*/", ""), ".gtc$", "") + ".vcf.gz",
      input_gtc = AutoCall.gtc_file,
      gender_gtc = GenderAutocall.gtc_file,
      extended_chip_manifest_file = extended_chip_manifest_file,
      cluster_file = cluster_file,
      normalization_manifest_file = AutoCall.bead_pool_manifest_csv_file,
      zcall_thresholds_file = zcall_thresholds_file,
      sample_alias = sample_alias,
      analysis_version_number = analysis_version_number,
      reported_gender = reported_gender,
      fingerprint_gender_file = fingerprint_gender_file,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      disk_size = disk_size,
      preemptible_tries = preemptible_tries
  }

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

  call CollectArraysVariantCallingMetrics {
    input:
      input_vcf_file = MergePedIntoVcf.output_vcf,
      input_vcf_index_file = MergePedIntoVcf.output_vcf_index,
      dbSNP_vcf_file = dbSNP_vcf,
      dbSNP_vcf_index_file = dbSNP_vcf_index,
      call_rate_threshold = call_rate_threshold,
      output_metrics_basename = sub(sub(AutoCall.gtc_file, "gs://.*/", ""), ".gtc$", ""),
      disk_size = disk_size,
      preemptible_tries = preemptible_tries
  }

  call VcfToIntervalList {
    input:
      vcf_file = MergePedIntoVcf.output_vcf,
      interval_list_file_name = sub(sub(MergePedIntoVcf.output_vcf, "gs://.*/", ""), ".vcf.gz$", "") + ".interval_list",
      disk_size = disk_size,
      preemptible_tries = preemptible_tries
  }

  call CheckFingerprint {
    input:
      input_vcf_file = MergePedIntoVcf.output_vcf,
      input_vcf_index_file = MergePedIntoVcf.output_vcf_index,
      genotypes_vcf_file = fingerprint_genotypes_vcf_file,
      genotypes_vcf_index_file = fingerprint_genotypes_vcf_index_file,
      haplotype_database_file = haplotype_database_file,
      observed_sample_alias = sub(sub(AutoCall.gtc_file, "gs://.*/", ""), ".gtc$", ""),
      expected_sample_alias = sample_alias,
      output_metrics_basename = sub(sub(AutoCall.gtc_file, "gs://.*/", ""), ".gtc$", ""),
      disk_size = disk_size,
      preemptible_tries = preemptible_tries
  }

  call GenotypeConcordance {
    input:
      call_vcf_file = MergePedIntoVcf.output_vcf,
      call_vcf_index_file = MergePedIntoVcf.output_vcf_index,
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

  output {
    AutoCall.gtc_file
    MergePedIntoVcf.*
    CollectArraysVariantCallingMetrics.*
    CheckFingerprint.*
    GenotypeConcordance.*
  }
}
