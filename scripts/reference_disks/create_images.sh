#!/bin/bash

# uncomment/comment toggle verbose debugging
#set -x

set -euo pipefail

BUCKET="gcp-public-data--broad-references"

usage() {
    cat <<USAGE >&2
Usage $0 [-f] <input_tsv>

Arguments:
* <input_tsv> should contain three columns without a header:
    Column 1: Name of the object under gs://$BUCKET (unique)
    Column 2: Size of the file in bytes
    Column 3: Image name prefix or 'none'

* -f (optional):
    Force a run of the commands instead of only performing a dry run

Output:
* <image_name>.manifest.conf
    A HOCON file describing the image

Dependencies:
* \`sudo\`
* \`gcloud\`
* \`jq\`
USAGE
    exit 1
}

# Date part to add to image names
DATE_STR=$(date +"%Y-%m-%d")

# The block size to divide images, aka all images will be a multiple of 5 (metric) gigabytes
BLOCK_GBS=5
BLOCK_SIZE=$((BLOCK_GBS * 10 ** 9))

DRY_RUN=true
EXISTING_RESOURCES=false
MISSING_DEPENDENCIES=false

error_msg() {
    echo "$(tput setaf 1)ERROR:$(tput sgr0)" "$@" >&2
}

# Convert an image name prefix to an image name
get_image_name() {
    local image_name_prefix

    image_name_prefix="$1"
    # Not 100% sure of the origin of this naming scheme. For now, the naming scheme is copied from
    # https://github.com/broadinstitute/firecloud-develop/blob/6705a3d7cccc9f536dff78a8359b6d4b586e32ab/base-configs/cromwell/cromwell-reference-images.conf#L3
    echo "$image_name_prefix-public-$DATE_STR"
}

# Check for already existing resources
check_existing() {
    local image_name

    image_name="$1"

    if [[ -n "$(
        gcloud compute images list --filter "name=('$image_name')" --format "value(name)"
    )" ]]; then
        error_msg "Existing image. Use \`gcloud compute images delete $image_name\`"
        EXISTING_RESOURCES=true
    fi

    if mount | grep "/mnt/disks/$image_name" &>/dev/null; then
        error_msg "Mounted disk. Use \`sudo umount /mnt/disks/$image_name\`"
        EXISTING_RESOURCES=true
    fi

    if [[ -n "$(
        gcloud compute instances list --filter "disks.deviceName=('$image_name')" --format="value(disks[].deviceName)"
    )" ]]; then
        error_msg "Attached disk. Use \`gcloud compute instances detach-disk $(hostname) --disk $image_name\`"
        EXISTING_RESOURCES=true
    fi

    if [[ -n "$(
        gcloud compute disks list --filter "name=('$image_name')" --format "value(name)"
    )" ]]; then
        error_msg "Existing disk. Use \`gcloud compute disks delete $image_name\`"
        EXISTING_RESOURCES=true
    fi

    if [[ -e "$image_name.manifest.conf" ]]; then
        error_msg "Existing manifest. Use \`rm '$image_name.manifest.conf'\`"
        EXISTING_RESOURCES=true
    fi
}

# Convert the retrieved base64 value to decimal since that's what CromwellRefdiskManifestCreatorApp does, slowly
# https://github.com/broadinstitute/cromwell/blob/65/CromwellRefdiskManifestCreator/src/main/java/org/broadinstitute/manifestcreator/CromwellRefdiskManifestCreatorApp.java#L116
get_cloud_hash() {
    local cloud_object

    cloud_object="$1"

    # URL encode the object using jq. Easier than trying to figure out which python2/python3 is available.
    # And we're using jq again just below to retrieve the crc32c from the output.
    # https://stackoverflow.com/questions/37309551/how-to-urlencode-data-into-a-url-with-bash-or-curl#comments-37309630
    cloud_object_encoded="$(jq --raw-input --raw-output @uri <<<"$cloud_object")"

    # No need to waste IO & CPU. The crc32c is already in the cloud
    # https://cloud.google.com/storage/docs/json_api
    hash_base64=$(
        curl \
            --location --fail --silent --show-error \
            "https://content-storage.googleapis.com/storage/v1/b/$BUCKET/o/$cloud_object_encoded?fields=crc32c" |
            jq --raw-output .crc32c
    )

    # Decode base64 to raw bytes, use od to convert to hex, then use awk to remove all whitespace
    # https://stackoverflow.com/questions/20262869/why-does-1-in-awk-print-the-current-line#answer-20263611
    hash_hex=$(
        base64 --decode <<<"$hash_base64" | od -t x1 -An | awk 'BEGIN { OFS=""; ORS="" } { gsub(/ /, "", $0) } 1'
    )

    # Convert the hex to decimal
    # https://stackoverflow.com/questions/13280131/hexadecimal-to-decimal-in-shell-script#answer-13280173
    echo $((16#$hash_hex))
}

maybe_cat_manifest() {
    local image_manifest
    local image_manifest_content

    image_manifest="$1"
    image_manifest_content="$2"

    if [[ "$DRY_RUN" == "false" ]]; then
        echo "adding to manifest: $image_manifest"
        echo "$image_manifest_content"
        echo "$image_manifest_content" >>"$image_manifest"
    else
        echo "DRY RUN to manifest: $image_manifest"
        echo "$image_manifest_content"
    fi
}

# Write the beginning of the manifest
start_manifest() {
    local image_manifest
    local image_name
    local image_size
    local image_identifier

    image_manifest="$1"
    image_name="$2"
    image_size="$3"
    image_identifier="projects/$(gcloud config get-value project)/global/images/$image_name"

    maybe_cat_manifest "$image_manifest" "{
    imageIdentifier: \"$image_identifier\"
    diskSizeGb: $image_size
    files: ["
}

# Write the closing section of the manifest
end_manifest() {
    local image_manifest

    image_manifest="$1"

    maybe_cat_manifest "$image_manifest" "    ]
}"
}

# Add a cloud object to the manifest
add_to_manifest() {
    local image_manifest
    local cloud_object
    local cloud_hash

    image_manifest="$1"
    cloud_object="$2"
    cloud_hash=$(get_cloud_hash "$cloud_object")

    maybe_cat_manifest "$image_manifest" "        {
            path: \"$BUCKET/$cloud_object\"
            crc32c: $cloud_hash
        }"
}

maybe_run() {
    if [[ "$DRY_RUN" == "false" ]]; then
        echo "running:" "$@"
        "$@"
    else
        echo "DRY RUN:" "$@"
    fi
}

# Return an image size based on the content size rounded up to the next block size
calculate_image_size() {
    local content_size
    local block_count
    local image_size

    content_size="$1"

    block_count=$(((content_size + BLOCK_SIZE - 1) / BLOCK_SIZE))
    # Minimum size is 10GB
    image_size=$((block_count * BLOCK_GBS < 10 ? 10 : block_count * BLOCK_GBS))

    echo "$image_size"
}

# Initialize and mount a persistent disk based on the image name
create_disk() {
    local image_name
    local image_size

    image_name="$1"
    image_size="$2"

    maybe_run gcloud compute disks create "$image_name" --size "${image_size}GB"

    maybe_run gcloud compute instances attach-disk "$(hostname)" --disk "$image_name" --device-name "$image_name"

    # copied from https://cloud.google.com/compute/docs/disks/add-persistent-disk#formatting
    maybe_run sudo mkfs.ext4 -m 0 -E lazy_itable_init=0,lazy_journal_init=0,discard "/dev/disk/by-id/google-$image_name"

    # copied from https://cloud.google.com/compute/docs/disks/add-persistent-disk#mounting
    maybe_run sudo mkdir -p "/mnt/disks/$image_name"
    maybe_run sudo mount -o discard,defaults "/dev/disk/by-id/google-$image_name" "/mnt/disks/$image_name"
}

# Download a file to the mounted persistent disk
download_file() {
    local image_name
    local cloud_path
    local local_path

    image_name="$1"
    cloud_path="gs://$BUCKET/$2"
    local_path="${cloud_path#"gs://"}"

    maybe_run sudo gsutil -m cp "$cloud_path" "/mnt/disks/$image_name/$local_path"
}

# Unmount the populated disk, copy the disk to a public image, then delete the disk
convert_disk_to_image() {
    local image_name

    image_name="$1"

    maybe_run sudo sync

    maybe_run sudo umount "/mnt/disks/$image_name"

    maybe_run sudo rmdir "/mnt/disks/$image_name"

    maybe_run gcloud compute instances detach-disk "$(hostname)" --disk "$image_name"

    maybe_run gcloud compute images create "$image_name" --source-disk "$image_name"

    # If you get an 403 error here then someone has shut off compute.images.setIamPolicy for the default compute SA.
    # Its then best to find out why they reverted those permissions and then consider running this script as another SA.
    # https://console.cloud.google.com/iam-admin/roles/details/projects%3Cbroad-dsde-cromwell-dev%3Croles%3CCustomComputeImagesAdmin?project=broad-dsde-cromwell-dev
    maybe_run gcloud compute images add-iam-policy-binding "$image_name" \
        --member "allAuthenticatedUsers" --role "roles/compute.imageUser"

    maybe_run gcloud --quiet compute disks delete "$image_name"
}

# Check that a dependent program is available
check_program() {
    local program

    program="$1"
    if ! command -v "$program" &>/dev/null; then
        error_msg "Program not found: \`$program\`"
        MISSING_DEPENDENCIES=true
    fi
}

# Check that a dependent file is readable
check_file() {
    local file

    file="$1"
    if [[ ! -r "$file" ]]; then
        error_msg "File not found: '$file'"
        MISSING_DEPENDENCIES=true
    fi
}

create_images() {
    local image_tsv
    local image_name_prefixes
    local image_name_prefix
    local image_name
    local image_manifest
    local content_size
    local image_size
    local cloud_objects
    local cloud_object

    image_tsv="$1"

    # Set the default zone for future gcloud commands to this machine's zone
    CLOUDSDK_COMPUTE_ZONE=$(
        gcloud compute instances list --filter "name=('$(hostname)')" --format "value(zone)"
    )
    export CLOUDSDK_COMPUTE_ZONE

    # Get a list of image name prefixes, excluding 'none'
    image_name_prefixes=$(cut -f 3 "$image_tsv" | grep -v '^none$' | sort -u)

    echo "Checking for existing resources..."
    for image_name_prefix in $image_name_prefixes; do
        image_name=$(get_image_name "$image_name_prefix")

        check_existing "$image_name"
    done

    if [[ "$EXISTING_RESOURCES" != "false" ]]; then
        exit 2
    fi

    echo "Creating images..."
    for image_name_prefix in $image_name_prefixes; do
        image_name=$(get_image_name "$image_name_prefix")
        image_manifest="$image_name.manifest.conf"

        # Get list of objects for this image
        cloud_objects=$(awk '$3 == "'"$image_name_prefix"'" { print $1 }' "$image_tsv")

        # Calculate the sum of the file contents
        content_size=$(
            # Use printf to handle large integers
            # https://stackoverflow.com/questions/450799/shell-command-to-sum-integers-one-per-line#answer-450821
            awk \
                '$3 == "'"$image_name_prefix"'" { gsub(/,/, "", $2); sum += $2 } END { printf "%.0f", sum }' \
                "$image_tsv"
        )

        # Calculate the image size
        image_size="$(calculate_image_size "$content_size")"

        # Create manifest header
        start_manifest "$image_manifest" "$image_name" "$image_size"
        # Add each item to the manifest
        for cloud_object in $cloud_objects; do
            add_to_manifest "$image_manifest" "$cloud_object"
        done
        # End the manifest
        end_manifest "$image_manifest"

        # Create and mount the persistent disk
        create_disk "$image_name" "$image_size"
        # Download files for each of the image names
        for cloud_object in $cloud_objects; do
            download_file "$image_name" "$cloud_object"
        done
        # Convert each of the disks to images and delete the disk
        convert_disk_to_image "$image_name"
    done

    echo "Done!"
}

if [[ "${1-}" == "-f" ]]; then
    DRY_RUN=false
    shift 1
fi

if [ "$#" -ne 1 ]; then
    usage
fi

IMAGE_TSV="$1"

check_program sudo
check_program gcloud
check_program jq
check_file "$IMAGE_TSV"

if [[ "$MISSING_DEPENDENCIES" != "false" ]]; then
    usage
fi

create_images "$IMAGE_TSV"
