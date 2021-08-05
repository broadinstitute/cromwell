version 1.0

task hello {
    input {
        Array[File] inputs
    }
    command {
        echo "Hello world!"
    }
    output {
        String out = read_string(stdout())
    }
    runtime {
      docker: "ubuntu:latest"
    }
}

task write_fofn {
    input {
        Int shard
    }

    Int num_inputs = 1700

    command <<<
        python <<CODE
        bucket = 'gs://bt-343-testing/'
        prefix = bucket + \
            'lorem-ipsum-dolor-sit-amet/consectetur-adipiscing-elit/nullam-in-aliquet-sapien/phasellus-at-feugiat-diam'

        output = open('inputs.txt', 'w')

        # Write an array of input files.
        # To replicate the scenario in BT-343 this should produce ~(6 * 285) or ~1700 inputs per shard.
        # Cycle through all of the inputs so a hash is requested for all of them (avoid the root workflow file hash
        # cache actor coalescing hash requests).

        lines = []
        for a in range(20):
            for b in range(10):
                for c in range(10):
                    lines.append(f'{prefix}/{a}-{"a"*64}/{b}-{"b"*64}/{c}-{"c"*64}/input.txt')

        x = (~{num_inputs} * ~{shard}) % ~{num_inputs}
        y = (~{num_inputs} * (~{shard} + 1)) % ~{num_inputs}

        if x < y:
            raw = lines[x:y]
        else:
            raw = lines[:y] + lines[x:]

        output.write('\n'.join(raw))

        output.close()
        CODE
    >>>
    output {
        Array[String] inputs = read_lines("inputs.txt")
    }
    runtime {
        docker: "python:latest"
    }
}

workflow lots_of_inputs_scattered {

    Int scatter_width = 700

    scatter (i in range(scatter_width)) {
        call write_fofn { input: shard = i }
    }

    scatter (i in range(scatter_width)) {
        call hello { input: inputs = write_fofn.inputs[i] }
    }

    output {
    }
}
