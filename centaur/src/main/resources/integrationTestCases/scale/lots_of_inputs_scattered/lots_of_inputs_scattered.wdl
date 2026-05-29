version 1.0

# Inspired by runs of the GATK SV WDL with ~700 shards localizing ~1700 inputs each.
# The inputs also had some specific characteristics that stressed Cromwell:
#
# * Long object paths of > 300 characters each, coming from subworkflows.
# * Largely unique, not taking much advantage of the RootWorkflowFileHashCacheActor and requiring
#   large numbers of hash requests to go out to Google.
# * Largely unique object path prefixes resulting in a large number of stanzas in the
#   `gcs_localization.sh` script, inflating the size of these scripts to ~1.7 MiB per shard.

task massive_localize {
    meta {
        volatile: false
    }
    input {
        Array[File] inputs
        String machine_type
    }
    command {
        echo "Hello world!"
    }
    output {
        String out = read_string(stdout())
    }
    runtime {
      docker: "ubuntu:latest"
      predefinedMachineType: machine_type
      maxRetries: 1
    }
}

task write_fofn {
    input {
        Int shard_index
        Int num_inputs
        String machine_type
    }

    command <<<
        python <<CODE
        bucket = 'gs://bt-343-testing/'
        prefix = bucket + \
            'lorem-ipsum-dolor-sit-amet/consectetur-adipiscing-elit/nullam-in-aliquet-sapien/phasellus-at-feugiat-diam'

        output = open('inputs_fofn.txt', 'w')

        # Write an array of input files.
        # To replicate the scenario in BT-343 this should produce ~(6 * 285) or ~1700 inputs per shard.
        # Cycle through all of the inputs so a hash is requested for all of them (avoid the root workflow file hash
        # cache actor coalescing hash requests).

        # Generate pool of 2,000 paths
        lines = []
        for a in range(20):
            for b in range(10):
                for c in range(10):
                    lines.append(f'{prefix}/{a}-{"a"*64}/{b}-{"b"*64}/{c}-{"c"*64}/input.txt')

        # Stagger each shard's window across the pool so different shards localize a different subset of inputs.
        pool_size = len(lines)
        start = (~{num_inputs} * ~{shard_index}) % pool_size
        end = start + ~{num_inputs}

        # Select within this pool (T) or wrap around to the next pool (F)
        if end <= pool_size:
            raw = lines[start:end]
        else:
            raw = lines[:end - pool_size] + lines[start:]

        output.write('\n'.join(raw))

        output.close()
        CODE
    >>>
    output {
        Array[String] inputs = read_lines("inputs_fofn.txt")
    }
    runtime {
        docker: "python:latest"
        predefinedMachineType: machine_type
        maxRetries: 1
    }
}

workflow lots_of_inputs_scattered {

    # An "expensive" T2D machine completes the task 3x faster and 50% cheaper than the default.
    input {
        Int scatter_width = 700
        Int num_inputs = 1700
        String machine_type = "t2d-standard-1"
    }

    scatter (i in range(scatter_width)) {
        call write_fofn {
            input:
                shard_index = i,
                num_inputs = num_inputs,
                machine_type = machine_type
        }
    }

    scatter (i in range(scatter_width)) {
        call massive_localize {
            input:
                inputs = write_fofn.inputs[i],
                machine_type = machine_type
        }
    }

    output {
    }
}
