workflow check_cromwell_root {
    call check_path
}

task check_path {

	runtime {
		docker: "ubuntu:latest"
	}

    command {
    	if [[ ! -e /cromwell_root ]]
    	then
      		echo "Could not find Cromwell root under /cromwell_root, exiting."
      		exit 1
    	fi
    }
}
