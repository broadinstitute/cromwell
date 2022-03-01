version 1.0

workflow bad_stdout {
    File f = stdout()

    output {
        String s = read_string(f)
    }
}
