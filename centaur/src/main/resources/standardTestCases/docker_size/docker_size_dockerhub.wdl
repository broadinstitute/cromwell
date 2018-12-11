version 1.0

task large_dockerhub_image {
    input {
        String docker_image
    }
    command {
        apt-get install --assume-yes jq > /dev/null
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
    call large_dockerhub_image as large_dockerhub_image_with_tag { input: docker_image = "broadinstitute/gatk:4.0.5.1" }
    # This is 4.0.5.2 - Use a different image to make sure we get the size for both
    call large_dockerhub_image as large_dockerhub_image_with_hash { input: docker_image = "broadinstitute/gatk@sha256:64c125136d6168250dfab99fd26215f40e6193a3cf1cd166fe70a4b8dd2d621a" }
}
