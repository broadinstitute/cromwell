#!/bin/bash
# The `papi_v2_log` Centaur test is opinionated about the number of log messages around localization/delocalization.
# The trace logging of `set -x` must be turned off for the `papi_v2_log` test to pass.
set +x
set -euo pipefail

gsutil_log=gsutil.log

NO_REQUESTER_PAYS_COMMAND=""
REQUESTER_PAYS_COMMAND=""


private::delocalize_file() {
  local cloud="$1"
  local container="$2"
  local rpflag="$3"
  local required="$4"
  local content="$5"

  # From Thibault:
  #
  # As per https://cloud.google.com/storage/docs/gsutil/addlhelp/HowSubdirectoriesWork, rule #2
  # If one attempts a
  #  gsutil cp /local/file.txt gs://bucket/subdir/file.txt
  #  AND
  #  there exists a folder gs://bucket/subdir/file.txt_thisCouldBeAnything
  #  then gs://bucket/subdir/file.txt will be treated as a directory, and /local/file.txt will be copied under gs://bucket/subdir/file.txt/file.txt
  #  and not gs://bucket/subdir/file.txt.
  #
  # By instead using the parent directory (and ensuring it ends with a slash), gsutil will treat that as a directory and put the file under it.
  # So the final gsutil command will look something like gsutil cp /local/file.txt gs://bucket/subdir/
  local cloud_parent=$(dirname "$cloud")
  cloud_parent="${cloud_parent}/"

  if [[ -f "$container" && -n "$content" ]]; then
    rm -f "$HOME/.config/gcloud/gce" && gsutil ${rpflag} -m -h "Content-Type: $content" cp "$container" "$cloud_parent" > "$gsutil_log" 2>&1
  elif [[ -f "$container" ]]; then
    rm -f "$HOME/.config/gcloud/gce" && gsutil ${rpflag} -m cp "$container" "$cloud_parent" > "$gsutil_log" 2>&1
  elif [[ -e "$container" ]]; then
    echo "File output '$container' exists but is not a file"
    exit 1
  elif [[ "$required" = "required" ]]; then
    echo "Required file output '$container' does not exist."
    exit 1
  fi
}


private::delocalize_directory() {
  local cloud="$1"
  local container="$2"
  local rpflag="$3"
  local required="$4"
  local content="$5"

  if [[ -d "$container" && -n "$content" ]]; then
    rm -f "$HOME/.config/gcloud/gce" && gsutil ${rpflag} -m -h "Content-Type: $content" rsync -r "$container" "$cloud" > "$gsutil_log" 2>&1
  elif [[ -d "$container" ]]; then
    rm -f "$HOME/.config/gcloud/gce" && gsutil ${rpflag} -m rsync -r "$container" "$cloud" > "$gsutil_log" 2>&1
  elif [[ -e "$container" ]]; then
    echo "Directory output '$container' exists but is not a directory"
    exit 1
  elif [[ "$required" = "required" ]]; then
    echo "Required directory output '$container' does not exist."
    exit 1
  fi
}


private::delocalize_file_or_directory() {
  local cloud="$1"
  local container="$2"
  local rpflag="$3"
  local required="$4"
  local content="$5"

  # required must be optional for 'file_or_directory' and was checked in the caller
  if [[ -f "$container" ]]; then
    private::delocalize_file "$cloud" "$container" "$rpflag" "$required" "$content"
  elif [[ -d "$container" ]]; then
    private::delocalize_directory "$cloud" "$container" "$rpflag" "$required" "$content"
  elif [[ -e "$container" ]]; then
    echo "'file_or_directory' output '$container' exists but is neither a file nor a directory"
    exit 1
  fi
}


private::timestamped_message() {
  printf '%s %s\n' "$(date -u '+%Y/%m/%d %H:%M:%S')" "$1"
}


private::localize_message() {
  local cloud="$1"
  local container="$2"
  local message=$(printf "Localizing input %s -> %s" "$cloud" "$container")
  private::timestamped_message "${message}"
}


private::delocalize_message() {
  local cloud="$1"
  local container="$2"
  local message=$(printf "Delocalizing output %s -> %s" "$container" "$cloud")
  private::timestamped_message "${message}"
}


# Requires both NO_REQUESTER_PAYS_COMMAND and USE_REQUESTER_PAYS_COMMAND to be set.
private::determine_requester_pays() {
  local max_attempts="$1"
  local attempt=1
  shift

  local command="$NO_REQUESTER_PAYS_COMMAND"
  local use_requester_pays=false
  # assume the worst
  USE_REQUESTER_PAYS=error

  while [[ ${attempt} -le ${max_attempts} ]]; do
    if eval ${command} > ${gsutil_log} 2>&1 ; then
      USE_REQUESTER_PAYS=${use_requester_pays}
      break
    elif [[ "$use_requester_pays" = "false" ]]; then
      if grep -q "Bucket is requester pays bucket but no user project provided." ${gsutil_log}; then
        use_requester_pays=true
        command="$REQUESTER_PAYS_COMMAND"
      else
        attempt=$((attempt + 1))
      fi
    else
      attempt=$((attempt + 1))
    fi
  done

  if [[ ${attempt} -gt ${max_attempts} ]]; then
    echo "Error attempting to localize file with command: '$command'"
    cat ${gsutil_log}
  fi
}


localize_files() {
  local project="$1"
  local max_attempts="$2"
  local container_parent="$3"
  local first_cloud_file="$4"
  shift 4

  local num_cpus=$(grep -c ^processor /proc/cpuinfo)
  # 32 is the max component count currently supported by gsutil cp.
  if [[ ${num_cpus} -gt 32 ]]; then
    num_cpus=32
  fi

  # We need to determine requester pays status of the first file attempting at most `max_attempts` times.
  NO_REQUESTER_PAYS_COMMAND="mkdir -p '$container_parent' && gsutil -o 'GSUtil:parallel_thread_count=1' -o 'GSUtil:sliced_object_download_max_components=${num_cpus}' cp '$first_cloud_file' '$container_parent'"
  REQUESTER_PAYS_COMMAND="gsutil -o 'GSUtil:parallel_thread_count=1' -o 'GSUtil:sliced_object_download_max_components=${num_cpus}' -u $project cp '$first_cloud_file' '$container_parent'"

  basefile=$(basename "$first_cloud_file")
  private::localize_message "$first_cloud_file" "${container_parent}${basefile}"
  private::determine_requester_pays ${max_attempts}

  if [[ ${USE_REQUESTER_PAYS} = true ]]; then
    rpflag="-u $project"
  elif [[ ${USE_REQUESTER_PAYS} = false ]]; then
    rpflag=""
  else
    # error
    exit 1
  fi


  if [[ $# -gt 0 ]]; then
    touch files_to_localize.txt
    while [[ $# -gt 0 ]]; do
      cloud="$1"
      basefile=$(basename "$cloud")
      container="${container_parent}${basefile}"
      private::localize_message "$cloud" "$container"
      echo "$cloud" >> files_to_localize.txt
      shift
    done

    attempt=1
    while [[ ${attempt} -le ${max_attempts} ]]; do
      # parallel transfer the remaining files
      if cat files_to_localize.txt | gsutil -o "GSUtil:parallel_thread_count=1" -o "GSUtil:sliced_object_download_max_components=${num_cpus}" -m ${rpflag} cp -I "$container_parent"; then
        break
      else
        attempt=$((attempt + 1))
      fi
    done
  fi
}


# Requires known requester pays status.
private::localize_directory() {
  local cloud="$1"
  local container="$2"
  local max_attempts="$3"
  local rpflag="$4"

  local attempt=1
  private::localize_message "$cloud" "$container"
  while [[ ${attempt} -lt ${max_attempts} ]]; do
    # Do not quote rpflag, when that is set it will be -u project which should be two distinct arguments.
    if mkdir -p "${container}" && rm -f "$HOME/.config/gcloud/gce" && gsutil ${rpflag} -m rsync -r "${cloud}" "${container}" > /dev/null 2>&1; then
      break
    else
      attempt=$(($attempt + 1))
    fi
  done

  if [[ ${attempt} -gt ${max_attempts} ]]; then
    exit 1
  fi
}


# Called from the localization script with unknown requester pays status on the source bucket. This attempts to localize
# the first input directory without requester pays. If that fails with a requester pays error, this attempts again with
# the project flag required for requester pays. Both no-requester-pays and requester-pays attempts are retried up to
# max_attempts times. Once requester pays status is determined via the first directory the remaining files are localized
# with or without the project flag as appropriate.
localize_directories() {
  local project="$1"
  local max_attempts="$2"
  local cloud_directory="$3"
  local container_directory="$4"
  shift 4

  BASE_COMMAND="private::localize_directory '${cloud_directory}' '${container_directory}' '${max_attempts}'"
  NO_REQUESTER_PAYS_COMMAND="${BASE_COMMAND} ''"
  REQUESTER_PAYS_COMMAND="${BASE_COMMAND} '-u $project'"

  private::determine_requester_pays ${max_attempts}

  if [[ ${USE_REQUESTER_PAYS} = true ]]; then
    rpflag="-u $project"
  elif [[ ${USE_REQUESTER_PAYS} = false ]]; then
    rpflag=""
  else
    exit 1
  fi

  while [[ $# -gt 0 ]]; do
    cloud_directory="$1"
    container_directory="$2"
    shift 2
    private::localize_directory "$cloud_directory" "$container_directory" "$max_attempts" "$rpflag"
  done
}


# Handles all delocalizations for a transfer bundle (a grouping of file, directories, or files_or_directories targeting
# a single GCS bucket).
delocalize() {
  local project="$1"
  local max_attempts="$2"

  shift 2

  # Whether the requester pays status of the GCS bucket is certain. rp status is presumed false until proven otherwise.
  local rp_status_certain=false
  local use_requester_pays=false

  while [[ $# -gt 0 ]]; do
    file_or_directory="$1"
    cloud="$2"
    container="$3"
    required="$4"
    content_type="$5"

    shift 5

    if [[ "$file_or_directory" = "file" ]]; then
      transfer_fn_name="private::delocalize_file"
    elif [[ "$file_or_directory" = "directory" ]]; then
      transfer_fn_name="private::delocalize_directory"
    elif [[ "$file_or_directory" = "file_or_directory" ]]; then
      transfer_fn_name="private::delocalize_file_or_directory"
    else
      echo "file_or_directory must be 'file' or 'directory' or 'file_or_directory' but got '$file_or_directory'"
      exit 1
    fi

    if [[ "$required" != "required" && "$required" != "optional" ]]; then
      echo "'required' must be 'required' or 'optional' but got '$required'"
      exit 1
    elif [[ "$required" = "required" && "$file_or_directory" = "file_or_directory" ]]; then
      echo "Invalid combination of required = required and file_or_directory = file_or_directory, file_or_directory only valid with optional secondary outputs"
      exit 1
    fi

    # Log what is being delocalized (at least one test depends on this).
    private::delocalize_message "$cloud" "$container"

    attempt=1
    transfer_rc=0
    # Loop attempting transfers for this file or directory while attempts are not exhausted.
    while [[ ${attempt} -le ${max_attempts} ]]; do

      if [[ ${use_requester_pays} = true ]]; then
        rpflag="-u ${project}"
      else
        rpflag=""
      fi

      # Note the localization versions of transfer functions are passed "required" and "content_type" parameters they will not use.
      if ${transfer_fn_name} "$cloud" "$container" "$rpflag" "$required" "$content_type"; then
        if [[ "$required" = "false" && ! -e "$container" ]]; then
          # Do not set rp_status_certain=true if an optional file was absent and no transfer was attempted.
          break
        else
          rp_status_certain=true
          break
        fi
      else
        private::timestamped_message "${transfer_fn_name} \"$cloud\" \"$container\" \"$rpflag\" \"$required\" \"$content_type\" failed"

        # Print the reason of the failure.
        cat "${gsutil_log}"

        # If the requester pays status of the GCS bucket is not certain look for requester pays errors.
        if [[ ${rp_status_certain} = false ]]; then
          if grep -q "Bucket is requester pays bucket but no user project provided." "${gsutil_log}"; then
            private::timestamped_message "Retrying with user project"
            use_requester_pays=true
            # Do not increment the attempt number, a requester pays failure does not count against retries.
            # Do mark that the bucket in question is certain to have requester pays status.
            rp_status_certain=true
          else
            # Requester pays status is not certain but this transfer failed for non-requester pays reasons.
            # Increment the attempt number.
            attempt=$((attempt+1))
          fi
        else
          attempt=$((attempt+1))
        fi
      fi
    done
    if [[ ${attempt} -gt ${max_attempts} ]]; then # out of attempts
      exit ${transfer_rc}
    fi
  done

  rm -f "${gsutil_log}"
}


# Required for files whose names are not consistent between cloud and container. There should be very few of these,
# the monitoring script being the single known example.
localize_singleton_file() {
  local project="$1"
  local max_attempts="$2"
  local cloud="$3"
  local container="$4"

  local container_parent=$(dirname "$container")

  private::localize_message "$cloud" "$container"
  NO_REQUESTER_PAYS_COMMAND="mkdir -p '$container_parent' && gsutil cp '$cloud' '$container'"
  REQUESTER_PAYS_COMMAND="gsutil -u $project cp '$cloud' '$container'"
  # As a side effect of determining requester pays this one file will be localized.
  private::determine_requester_pays ${max_attempts}
}
