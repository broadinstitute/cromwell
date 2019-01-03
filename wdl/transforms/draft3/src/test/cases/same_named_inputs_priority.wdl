version 1.0

# Check that we do not "leak" references to `out` between scopes - i.e. should never get "hello1234"
workflow same_named_inputs_priority {
    String out = "hello"
    call echo as a {
        input:
            out = out + "1"
    }
    call echo as b {
        input:
            out = out + "2"
    }
    call echo as c {
        input:
            out = out + "3"
    }
    call echo as d {
        input:
            out = out + "4"
    }
}
task echo {
    runtime {
      docker: "ubuntu:latest"
    }
    input {
        String out
    }
    command {
        echo ~{out}
    }
    output {
        String result = out
    }
}
