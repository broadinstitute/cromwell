version 1.0

task ls {
    input {
        Array[File] files
    }
    command {
        ls ${sep=" " files}
    }

    output {
        Array[String] filename = read_lines(stdout())
    }
}


workflow array_autoboxing_inputs {
    input {
        Array[File] files
    }

    scatter (i in files) {
        call ls {
            input: files = i
        }
    }

    output {
        Array[Array[String]] F = ls.filename
    }
}
