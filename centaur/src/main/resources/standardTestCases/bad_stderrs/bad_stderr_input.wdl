version 1.0

workflow bad_stderr {

    input {
        File f = stderr()
    }
    output {
        String s = read_string(f)
    }
}
