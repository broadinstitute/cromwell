version 1.0

workflow DownloadFromHA {

    String data_url = "gs://bucket/random_file.txt"

    call UseBadBasename { input: data_url = data_url }
}

task UseBadBasename {
    input {
        String data_url
    }

    String file_path = basename(sub(data_url, "?.*", ""))

    command <<<
        echo ~{file_path}
    >>>

    output {
        File out = stdout()
    }

    runtime {
        docker: "ubuntu:latest"
    }
}
