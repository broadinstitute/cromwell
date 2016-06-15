task centaur {
    String branch
    File conf
    File pem

    command<<<
        sudo apt-get update
        sudo apt-get --yes install git
        git clone https://github.com/broadinstitute/centaur.git
        centaur/test_cromwell.sh -b${branch} -c${conf}
    >>>

    output {
        String centaur_log = read_string(stdout())
    }

    runtime {
        docker: "ubuntu:latest"
        # FIXME: Some bigass machine here
    }
}

workflow centaur {
    call centaur
}