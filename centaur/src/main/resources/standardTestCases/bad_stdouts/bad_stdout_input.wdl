version 1.0

workflow bad_stdout {
    input {
        File f = stdout()
    }
    output {
        String s = read_string(f)
    }
}
