version 1.0

workflow bad_stdout {
    output {
        File f = stdout()
    }
}
