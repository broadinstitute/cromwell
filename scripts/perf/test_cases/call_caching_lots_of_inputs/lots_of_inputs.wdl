version 1.0

task lots_of_inputs {
    input {
        Array[File] files
    }
    
    command {
        echo ~{sep=' ' files}
    }
    
    runtime {
        docker: "ubuntu@sha256:de774a3145f7ca4f0bd144c7d4ffb2931e06634f11529653b23eba85aef8e378"
    }
    
    output {
        File out = stdout()
    }
}

workflow call_cache_read_input_files_perf_test {
    input {
        Array[File] files
    }
    scatter(i in range(100)) {
        call lots_of_inputs { input: files = files }
    }
}
