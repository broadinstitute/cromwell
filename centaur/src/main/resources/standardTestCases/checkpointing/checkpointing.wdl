version 1.0

workflow checkpointing {
  call count { input: count_to = 100 }
  output {
    String preempted = count.preempted
  }
}

task count {
  input {
    Int count_to
  }

  meta {
    volatile: true
  }

  command <<<
    # Read from the my_checkpoint file if there's content there:
    FROM_CKPT=$(cat my_checkpoint | tail -n1 | awk '{ print $1 }')
    FROM_CKPT=${FROM_CKPT:-1}

    # We don't want any single VM run the entire count, so work out the max counter value for this attempt:
    MAX="$(($FROM_CKPT + 66))"

    INSTANCE_NAME=$(curl -s "http://metadata.google.internal/computeMetadata/v1/instance/name" -H "Metadata-Flavor: Google")
    echo "Discovered instance: $INSTANCE_NAME"

    # Run the counter:
    echo '--' >> my_checkpoint
    for i in $(seq $FROM_CKPT ~{count_to})
    do
      echo $i
      echo $i ${INSTANCE_NAME} $(date) >> my_checkpoint

      # If we're over our max, simulate "preempting" the VM by killing it:
      if [ "${i}" -gt "${MAX}" ]
      then
        fully_qualified_zone=$(curl -s -H "Metadata-Flavor: Google" http://metadata.google.internal/computeMetadata/v1/instance/zone)
        zone=$(basename "$fully_qualified_zone")
        gcloud compute instances delete ${INSTANCE_NAME} --zone=$zone -q
      fi

      sleep 1
    done

    # Prove that we got preempted at least once:
    FIRST_INSTANCE=$(cat my_checkpoint | head -n1 | awk '{ print $2 }')
    LAST_INSTANCE=$(cat my_checkpoint | tail -n1 | awk '{ print $2 }')
    if [ "${FIRST_INSTANCE}" != "LAST_INSTANCE" ]
    then
      echo "GOTPREEMPTED" > preempted.txt
    else
      echo "NEVERPREEMPTED" > preempted.txt
    fi
  >>>

  runtime {
    docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:slim"
    preemptible: 3
    checkpointFile: "my_checkpoint"
  }

  output {
    File checkpoint_log = "my_checkpoint"
    String preempted = read_string("preempted.txt")
  }
}
