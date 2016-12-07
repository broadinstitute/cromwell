task centaur {
    String centaur_branch
    File conf
    File pem
    File cromwell_jar
    File token
    String secret = read_string(token)

    command<<<
        mkdir -p /cromwell_root/tmp/ivy2
        export SBT_OPTS=-Dsbt.ivy.home=/cromwell_root/tmp/.ivy2
        git clone https://github.com/broadinstitute/centaur.git
        cd centaur
        git checkout ${centaur_branch}
        cd ..
        centaur/test_cromwell.sh -j${cromwell_jar} -c${conf} -r/cromwell_root -t ${secret} -elocaldockertest -p100
    >>>

    output {
       File cromwell_log = "/cromwell_root/logs/cromwell.log"
       File centaur_log_file = "/cromwell_root/logs/centaur.log"
    }

    runtime {
        docker: "geoffjentry/centaur-cromwell:latest" 
        cpu: "8"
        zones: "us-central1-b"
        failOnStderr: false
    }
}
workflow centaur_workflow {
    call centaur
}
