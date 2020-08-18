version 1.0

task large_dockerhub_image {
    input {
        String docker_image
    }
    command {
        which jq > /dev/null || (apt-get update > /dev/null && apt-get install -y jq > /dev/null)
        NAME=`curl -s -H "Metadata-Flavor: Google" http://metadata.google.internal/computeMetadata/v1/instance/name`
        ZONE=`basename \`curl -s -H "Metadata-Flavor: Google" http://metadata.google.internal/computeMetadata/v1/instance/zone\``
        PROJECT=`curl -s -H "Metadata-Flavor: Google" http://metadata.google.internal/computeMetadata/v1/project/project-id`
        curl -s -H "Authorization: Bearer `gcloud auth print-access-token`" "https://www.googleapis.com/compute/v1/projects/$PROJECT/zones/$ZONE/disks/$NAME?fields=sizeGb" | jq -r '.sizeGb'
    }

    runtime {
        docker: docker_image
    }
    
    output {
        Int bootDiskSize = read_int(stdout())
    }
}

workflow docker_size_dockerhub {
    call large_dockerhub_image as large_dockerhub_image_with_tag { input: docker_image = "broadinstitute/gatk:4.1.8.0" }
    # This is 4.1.8.1 - Use a different image to make sure we get the size for both
    call large_dockerhub_image as large_dockerhub_image_with_hash { input: docker_image = "broadinstitute/gatk@sha256:8051adab0ff725e7e9c2af5997680346f3c3799b2df3785dd51d4abdd3da747b" }
}
