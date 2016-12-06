task url_grab {
    String url
    command <<<
        curl -# ${url} > news
        sleep 2
    >>>
    output {
        File news = "news"
        File progressBar = stderr()
    }
    runtime {
        docker: "tutum/curl:latest"
        failOnStderr: false
    }
}

task news_reader {
    File news

    command { 
      wc -l < ${news}
      sleep 2
    }
    output { Int news_size = read_int(stdout()) }
    runtime { docker: "ubuntu:latest" }
}

workflow curl_wf {
    call url_grab as newsgrab { input: url = "http://www.bbc.com/news" }
    call news_reader { input: news = newsgrab.news }
}
