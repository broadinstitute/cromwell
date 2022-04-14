version 1.0

workflow bad_stderr {
    output {
        File f = stderr()
    }
}
