version 1.0


workflow ReferenceDiskMountedOnlyIfRequested {
    call MentionsNirvanaReference {}
    output {
        Boolean disk_mounted = MentionsNirvanaReference.disk_mounted
    }
}


task MentionsNirvanaReference {
    input {
        # Tiny 116 byte reference file, certainly not worth attaching a 55 GiB Nirvana reference disk for this.
        File mention =
            "gs://broad-public-datasets/gvs/vat-annotations/Nirvana/3.18.1/SupplementaryAnnotation/GRCh38/MITOMAP_20200819.nsa.idx"
    }
    command <<<
        PS4='\D{+%F %T} \w $ '
        set -o nounset -o xtrace

        # Debug output
        lsblk > lsblk.out

        # Look for mounted disks other than the main /mnt/disks/cromwell_root volume
        CANDIDATE_MOUNT_POINT=$(lsblk | grep -E -v '/mnt/disks/cromwell_root\>' | sed -E -n 's!.*(/mnt/disks/[^/]+).*!\1!p')

        if [[ ! -z ${CANDIDATE_MOUNT_POINT} ]]; then
            echo "Found unexpected mounted disk, investigating further."
            find ${CANDIDATE_MOUNT_POINT} -print | tee find.out

            if grep -i nirvana find.out; then
                echo "Found what appears to be a Nirvana reference disk."
            else
                echo "Found unknown volume mounted, see 'find.out' for manifest."
            fi
            echo true > disk_mounted.out
        else
            echo false > disk_mounted.out
        fi
    >>>
    runtime {
        docker: "ubuntu:latest"
        backend: "GCPBATCH-Reference-Disk-Localization"
    }
    output {
        File lsblk = "lsblk.out"
        Boolean disk_mounted = read_boolean("disk_mounted.out")
        File? find_out = "find.out"
    }
}
