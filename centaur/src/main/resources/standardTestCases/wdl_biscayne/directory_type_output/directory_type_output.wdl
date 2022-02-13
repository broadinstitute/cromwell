version development

# CROM-6875 repro WDL to exercise Directory outputs courtesy Giulio Genovese. Since the Directory type does not exist in
# WDL versions 1.0 or draft-2, the bug this is checking for cannot and does not exist in those WDL versions.
workflow main {
    call main { input: s1 = "x", s2 = "y" }
    output { Array[File] f = main.f }
}

task main {
    input {
        String s1
        String s2
    }

    command <<<
        set -euo pipefail
        mkdir d
        touch "d/~{s1}"
        touch "d/~{s2}"
        echo -e "d/~{s1}\nd/~{s2}"
    >>>

    output {
        Directory d = "d"
        Array[File] f = read_lines(stdout())
    }

    runtime {
        docker: "debian:stable-slim"
    }
}
