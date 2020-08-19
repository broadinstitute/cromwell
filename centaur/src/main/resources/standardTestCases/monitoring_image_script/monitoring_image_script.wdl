version 1.0

workflow monitoring_image_script {
    call create_file_via_monitoring_image_script

    output {
        String monitoring_script_created = create_file_via_monitoring_image_script.monitoring_script_created
    }
}

task create_file_via_monitoring_image_script {
    meta {
        volatile: true
    }
    command {
        sleep 10
    }
    output {
        String monitoring_script_created = read_string("monitoring_script_created")
    }
    runtime {
        docker: "ubuntu"
    }
}
