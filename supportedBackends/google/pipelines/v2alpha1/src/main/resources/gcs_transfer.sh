#!/bin/bash
# The `papi_v2_log` Centaur test is opinionated about the number of log messages around localization/delocalization.
# The trace logging of `set -x` must be turned off for the `papi_v2_log` test to pass.
set +x

gsutil_log="gsutil_output.txt"


localize_file() {
  local cloud="$1"
  local container="$2"
  local rpflag="$3"
  # Do not quote rpflag, when that is set it will be -u project which should be two distinct arguments.
  rm -f "$HOME/.config/gcloud/gce" && gsutil ${rpflag} -m cp "$cloud" "$container" > "$gsutil_log" 2>&1
}

localize_directory() {
  local cloud="$1"
  local container="$2"
  local rpflag="$3"
  # Do not quote rpflag, when that is set it will be -u project which should be two distinct arguments.
  mkdir -p "${container}" && rm -f "$HOME/.config/gcloud/gce" && gsutil ${rpflag} -m rsync -r "${cloud}" "${container}" > "$gsutil_log" 2>&1
}

delocalize_file() {
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

delocalize_directory() {
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

delocalize_file_or_directory() {
  local cloud="$1"
  local container="$2"
  local rpflag="$3"
  local required="$4"
  local content="$5"

  # required must be optional for 'file_or_directory' and was checked in the caller
  if [[ -f "$container" ]]; then
    delocalize_file "$cloud" "$container" "$rpflag" "$required" "$content"
  elif [[ -d "$container" ]]; then
    delocalize_directory "$cloud" "$container" "$rpflag" "$required" "$content"
  elif [[ -e "$container" ]]; then
    echo "'file_or_directory' output '$container' exists but is neither a file nor a directory"
    exit 1
  fi
}

timestamped_message() {
  printf '%s %s\n' "$(date -u '+%Y/%m/%d %H:%M:%S')" "$1"
}

localize_message() {
  local cloud="$1"
  local container="$2"
  local message=$(printf "Localizing input %s -> %s" "$cloud" "$container")
  timestamped_message "${message}"
}

delocalize_message() {
  local cloud="$1"
  local container="$2"
  local message=$(printf "Delocalizing output %s -> %s" "$container" "$cloud")
  timestamped_message "${message}"
}

# Transfer a bundle of files or directories to or from the same GCS bucket.
transfer() {
  local direction="$1"
  local project="$2"
  local max_attempts="$3"

  shift 3 # direction + project + max_attempts

  if [[ "$direction" != "localize" && "$direction" != "delocalize" ]]; then
    echo "direction must be 'localize' or 'delocalize' but got '$direction'"
    exit 1
  fi

  # Whether the requester pays status of the GCS bucket is certain. rp status is presumed false until proven otherwise.
  local rp_status_certain=false
  local use_requester_pays=false

  local message_fn="${direction}_message"

  # Loop while there are still items in the bundle to transfer.
  while [[ $# -gt 0 ]]; do
    file_or_directory="$1"
    cloud="$2"
    container="$3"

    if [[ "$file_or_directory" = "file" ]]; then
      transfer_fn_name="${direction}_file"
    elif [[ "$file_or_directory" = "directory" ]]; then
      transfer_fn_name="${direction}_directory"
    elif [[ "$direction" = "delocalize" && "$file_or_directory" = "file_or_directory" ]]; then
      transfer_fn_name="delocalize_file_or_directory"
    else
      echo "file_or_directory must be 'file' or 'directory' or (for delocalization only) 'file_or_directory' but got '$file_or_directory' with direction = '$direction'"
      exit 1
    fi

    content_type=""
    required=""
    if [[ "${direction}" = "delocalize" ]]; then
      # 'required' and 'content type' only appear in delocalization bundles.
      required="$4"
      content_type="$5"
      if [[ "$required" != "required" && "$required" != "optional" ]]; then
        echo "'required' must be 'required' or 'optional' but got '$required'"
        exit 1
      elif [[ "$required" = "required" && "$file_or_directory" = "file_or_directory" ]]; then
        echo "Invalid combination of required = required and file_or_directory = file_or_directory, file_or_directory only valid with optional secondary outputs"
        exit 1
      fi
      shift 2 # required + content_type
    fi
    shift 3 # file_or_directory + cloud + container

    # Log what is being localized or delocalized (at least one test depends on this).
    ${message_fn} "$cloud" "$container"

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
      ${transfer_fn_name} "$cloud" "$container" "$rpflag" "$required" "$content_type"
      transfer_rc=$?

      if [[ ${transfer_rc} = 0 ]]; then
        rp_status_certain=true
        break
      else
        timestamped_message "${transfer_fn_name} \"$cloud\" \"$container\" \"$rpflag\" \"$required\" \"$content_type\" failed"

        # Print the reason of the failure.
        cat "${gsutil_log}"

        # If the requester pays status of the GCS bucket is not certain look for requester pays errors.
        if [[ ${rp_status_certain} = false ]]; then
          if grep -q "Bucket is requester pays bucket but no user project provided." "${gsutil_log}"; then
            timestamped_message "Retrying with user project"
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
}
