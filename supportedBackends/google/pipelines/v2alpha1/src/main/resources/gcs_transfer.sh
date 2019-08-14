#!/bin/bash
# The `papi_v2_log` Centaur test is opinionated about the number of log messages around localization/delocalization.
# The trace logging of `set -x` must be turned off for the `papi_v2_log` test to pass.
set +x

gsutil_log=$(mktemp /tmp/gsutil.XXXXXXXXXXXXXXXX)


private::localize_file() {
  local cloud="$1"
  local container="$2"
  local rpflag="$3"
  # Do not quote rpflag, when that is set it will be -u project which should be two distinct arguments.
  rm -f "$HOME/.config/gcloud/gce" && gsutil ${rpflag} -m cp "$cloud" "$container" > "$gsutil_log" 2>&1
}

private::localize_directory() {
  local cloud="$1"
  local container="$2"
  local rpflag="$3"
  # Do not quote rpflag, when that is set it will be -u project which should be two distinct arguments.
  mkdir -p "${container}" && rm -f "$HOME/.config/gcloud/gce" && gsutil ${rpflag} -m rsync -r "${cloud}" "${container}" > "$gsutil_log" 2>&1
}

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

# Transfer a bundle of files or directories to or from the same GCS bucket.
transfer() {
  # Begin the transfer with uncertain requester pays status and first attempting transfers without requester pays.
  private::transfer false false "$@"
}

private::transfer() {
  local rp_status_certain="$1"
  local use_requester_pays="$2"
  local direction="$3"
  local project="$4"
  local max_attempts="$5"

  shift 5 # rp_status_certain + use_requester_pays + direction + project + max_attempts

  if [[ "$direction" != "localize" && "$direction" != "delocalize" ]]; then
    echo "direction must be 'localize' or 'delocalize' but got '$direction'"
    exit 1
  fi

  # Whether the requester pays status of the GCS bucket is certain. rp status is presumed false until proven otherwise.

  local message_fn="private::${direction}_message"

  # If requester pays status is unknown, loop through the items in the transfer bundle until requester pays status is determined.
  # Once determined, the remaining items can be transferred in parallel.
  while [[ $# -gt 0 ]]; do
    file_or_directory="$1"
    cloud="$2"
    container="$3"

    if [[ "$file_or_directory" = "file" ]]; then
      transfer_fn_name="private::${direction}_file"
    elif [[ "$file_or_directory" = "directory" ]]; then
      transfer_fn_name="private::${direction}_directory"
    elif [[ "$direction" = "delocalize" && "$file_or_directory" = "file_or_directory" ]]; then
      transfer_fn_name="private::delocalize_file_or_directory"
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

      # Do not set rp_status_certain=true if an optional file was absent and no transfer was attempted.
      if [[ ${transfer_rc} = 0 && "$required" = "false" && ! -e "$container" ]]; then
        break
      elif [[ ${transfer_rc} = 0 ]]; then
        rp_status_certain=true
        break
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
