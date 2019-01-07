version 1.0

task shouldCompleteFast {
    input {
      Int a
    }
    command {
        echo "The number was: ~{a}"
    }
    output {
        Int echo = a
    }
    runtime {
      docker: "ubuntu:latest"
    }
}

task shouldCompleteSlow {
    input {
      Int a
    }
    command {
        echo "The number was: ~{a}"
        # More than 1 so this should finish second
        sleep 2
    }
    output {
        Int echo = a
    }
    runtime {
      docker: "ubuntu:latest"
    }
}

task failMeSlowly {
    input {
      Int a
    }
    command {
        echo "The number was: ~{a}"
        # Less than 2 so this should finish first
        sleep 1
        ./NOOOOOO
    }
    output {
        Int echo = a
    }
    runtime {
      docker: "ubuntu:latest"
    }
}

task shouldNeverRun {
    input {
      Int a
      Int b
    }
    command {
        echo "You can't fight in here - this is the war room ~{a + b}"
    }
    output {
        Int echo = a
    }
    runtime {
      docker: "ubuntu:latest"
    }
}

workflow wf {
    call shouldCompleteFast as A { input: a = 5 }
    call shouldCompleteFast as B { input: a = 5 }

    call failMeSlowly as ohNOOOOOOOO { input: a = A.echo }
    call shouldCompleteSlow as C { input: a = B.echo }

    call shouldNeverRun as D { input: a = ohNOOOOOOOO.echo, b = C.echo }
    call shouldCompleteSlow as E { input: a = C.echo }
}
