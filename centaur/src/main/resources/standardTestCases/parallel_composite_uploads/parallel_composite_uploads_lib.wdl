version 1.0

# This task and the following differ only by the `backend` runtime attribute but unfortunately that must be a literal.
task make_a_fake_bam_parallel_composite_uploads_off_in_config {
    command {
      # 1 GiB fake BAM
      dd if=/dev/urandom of=fake.bam bs=$((2 ** 10)) count=$((2 ** 20))
    }
    output {
      File bam = "fake.bam"
    }
    runtime {
      docker: "ubuntu:latest"
      backend: "Papiv2"
    }
}

task make_a_fake_bam_parallel_composite_uploads_on_in_config {
    command {
      # 1 GiB fake BAM
      dd if=/dev/urandom of=fake.bam bs=$((2 ** 10)) count=$((2 ** 20))
    }
    output {
      File bam = "fake.bam"
    }
    runtime {
      docker: "ubuntu:latest"
      backend: "Papiv2ParallelCompositeUploads"
    }
}

task check_composite {
   input {
     String bam_path
   }
   command <<<
     # Look for a Component-Count from a `gsutil stat` of the previously uploaded "bam".
     # This property should not exist if the file was not parallel composite uploaded.
     count=$(gsutil stat ~{bam_path} | grep 'Component-Count' | sed -E 's![^[:digit:]]+(.*)!\1!')
     if [[ -z "$count" ]]; then echo false; elif [[ "$count" -gt 1 ]]; then echo true; else echo false; fi
   >>>
   output {
     Boolean composite = read_boolean(stdout())
   }
   runtime {
     docker: "google/cloud-sdk:latest"
   }
}

