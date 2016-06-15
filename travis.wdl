task centaur {
    String branch
    File conf
    File pem

    command<<<
        git clone https://github.com/broadinstitute/centaur.git
        git clone https://github.com/broadinstitute/cromwell.git
        cd cromwell
        echo "cd"
        git checkout ${branch}
        echo "checkout"
        git pull
        echo "pull"
        sbt clean
        echo "clean"
        sbt assembly
        echo "assembly"
        java -Dconfig.file=${conf} -jar target/scala-2.11/cromwell-*.jar
        echo "done"
        exit 0
    >>>

    output {
        String centaur_log = read_string(stdout())
    }

    runtime {
        docker: "geoffjentry/centaur-cromwell:latest"
        cpu: "8"
        memory: "10 GB"
    }
}

workflow centaur {
    call centaur
}