version 1.0

task hello_memory_unit {
    input {
        String memory_amount
    }

    command {
        echo "hello ~{memory_amount}"
    }

    output {
        String out = read_string(stdout())
    }

    runtime {
        docker: "ubuntu"
        memory: memory_amount
    }
}

workflow memory_units {
    Array[String] memory_amounts = [
        "1000000000 B",
        "1000000 KiB",
        "1000000 Ki",
        "1000000 KB",
        "1000000 K",
        "1000 MiB",
        "1000 Mi",
        "1000 MB",
        "1000 M",
        "1 GiB",
        "1 Gi",
        "1 GB",
        "1 G",
        "0.001 TiB",
        "0.001 Ti",
        "0.001 TB",
        "0.001 T",
     ]

    scatter (memory_amount in memory_amounts) {
        call hello_memory_unit {
            input: memory_amount = memory_amount
        }
    }

    output {
        Array[String] outs = hello_memory_unit.out
    }
}
