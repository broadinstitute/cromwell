version 1.0

task large_gcr_image {
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

workflow docker_size_gcr {
   call large_gcr_image as large_gcr_image_with_tag { input: docker_image = "gcr.io/broad-dsde-cromwell-dev/centaur/gatk:4.0.5.1" }
   # This is 4.0.5.2 - Use a different image to make sure we get the size for both
   call large_gcr_image as large_gcr_image_with_hash { input: docker_image = "gcr.io/broad-dsde-cromwell-dev/centaur/gatk@sha256:d78b14aa86b42638fe2844def82816d002a134cc19154a21dac7067ecb3c7e06" }
}
