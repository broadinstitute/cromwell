version 1.0

workflow input_starved {
    input {
        String feedme
        String please
    }
    call echo as echoFeedMe {input: string=feedme}
    call echo as echoPlease {input: string=please}

    output {}
}

task echo {
    input {String string}
    command <<<echo ~{string}>>>
    output {String out = stdout()}
    runtime {docker: "debian:buster-slim"}
}