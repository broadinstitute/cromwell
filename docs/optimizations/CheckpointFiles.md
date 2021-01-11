# The 'Checkpoint File' Optimization

## Overview

Available in Cromwell 55 and higher.

This optimization hopes to resolve the issue of your worker VM being preempted 9 hours and 55 minutes into the runtime of
a 10 hour job and having no option but to re-run the entire computation again.

### Description

Specifying a `checkpointFile` value in a task's `runtime` section designates a checkpoint file which will periodically be
copied to cloud storage every 10 minutes.
This checkpoint file will then be restored automatically on subsequent attempts if the job is interrupted.

**Note:** Although the checkpoint file is deleted if the task succeeds, additional charges may accrue storing the checkpoint file during
the running of the task, if the task is aborted or otherwise stopped externally, and by transferring it between the VM and the cloud. These
cost should be minor, especially balanced against the performance and cost benefits of being able to restore from the
checkpoint when preemptible VMs are interrupted.   

### Effect on Call Caching

The presence or absence of the `checkpointFile` attribute is not considered when determining whether to call cache.  

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
    # Note: Cromwell will stage the checkpoint file on recovery attempts.
    # This task checks the 'my_checkpoint' file for a counter value, or else
    # initializes the counter at '1':
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
    # Note: This checkpointFile attribute is what signals to Cromwell to save
    # the designated checkpoint file:
    checkpointFile: "my_checkpoint"
  }

  output {
    # Note: This task also uses the checkpoint as its output. This is not
    # required for checkpointing to work:
    Array[String] out = read_lines("my_checkpoint")
  }
}
```

## Backend Support

Cromwell supports the `checkpointFile` attribute on the following backends:

* The Google PAPIv2 (alpha1) backend
* The Google Life Sciences (beta) backend
