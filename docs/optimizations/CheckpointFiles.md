# The 'Checkpoint File' Optimization

## Overview

Available in Cromwell version 55 and higher.

### Description

Specifying a `checkpointFile` value in a task's `runtime` section designates a checkpoint file which will occasionally be
copied out to cloud storage. This file will then be restored automatically on subsequent attempts if the job is interrupted.

*Note:* Additional charges may accrue storing the checkpoint file, and by transferring it to and from the cloud. This
should be balanced against the performance and cost benefits of being able to restore from the checkpoint when preemptible VMs are interrupted.   

### Effect on Call Caching

The presence or absence of the `checkpointingFile` attribute is not considered when determining whether to call cache.  

### Example

The following WDL demonstrates the use of the `checkpointFile` optimization. It has a command which is checkpoint-aware:

* It starts by attempting to restore state from the `my_checkpoint` file (or starts at `1` if the checkpoint is empty)
* Then it counts up to 100, printing out the current counter value and a date timestamp at each value.

To make the checkpointing work, the `runtime` section specifies `checkpointFile: "my_checkpoint"`.

```wdl
version 1.0

workflow count_wf {
  call count { input: count_to = 100 }
}

task count {
  input {
    Int count_to
  }

  command <<<
    # Read from the my_checkpoint file if there's content there:
    FROM_CKPT=$(cat my_checkpoint | tail -n1 | awk '{ print $1 }')
    FROM_CKPT=${FROM_CKPT:-1}

    echo '--' >> my_checkpoint
    for i in $(seq $FROM_CKPT ~{count_to})
    do
      echo $i $(date) >> my_checkpoint
      sleep 4
    done
  >>>

  runtime {
    docker: "ubuntu:latest"
    preemptible: 3
    checkpointFile: "my_checkpoint"
  }

  output {
    Array[String] out = read_lines("my_checkpoint")
  }
}
```

## Backend Support

The `checkpointFile` attribute is currently understood by and implemented for:

* The Google PAPIv2 (alpha1) backend
* The Google Life Sciences (beta) backend
