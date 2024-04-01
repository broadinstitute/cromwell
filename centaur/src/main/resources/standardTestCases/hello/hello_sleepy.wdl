task takeanap {
  Int sleepytime
  command {
    echo "Time to sleep for ${sleepytime} seconds"
    sleep ${sleepytime}
    echo "Ahh, I feel so refreshed!"
  }
  output {
    String naplog = read_string(stdout())
  }
  runtime {
    docker: "ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950"
  }
}

workflow wf_hellosleepy {
  call takeanap
  output {
     takeanap.naplog
  }
}
