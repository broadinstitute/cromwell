task centaur {
    String cromwell_branch
    File conf
    File pem

    command<<<
        mkdir -p /cromwell_root/tmp/.ivy2
        export SBT_OPTS=-Dsbt.ivy.home=/cromwell_root/tmp/.ivy2
        cd /centaur
        git checkout develop
        git pull
        ./test_cromwell.sh -b${cromwell_branch} -r/cromwell_root -c../${conf} -p25
    >>>

    output {
        File cromwell_log = "/cromwell_root/logs/cromwell.log"
        File centaur_log_file = "/cromwell_root/logs/centaur.log"
    }

    runtime {
        docker: "geoffjentry/centaur-cromwell:latest"
        cpu: "32"
        memory: "10 GB"
        zones: "us-central1-b"
        failOnStderr: false
    }
}

workflow centaur {
    call centaur
}