workflow check_cromwell_root {
    call print_cromwell_root
}

task print_cromwell_root {

	runtime {
		docker: "ubuntu:latest"
	}

    command {
    	if [[ -e /cromwell_root/gcs_delocalization.sh ]]
    	then
      		CROMWELL_ROOT=/cromwell_root
    	elif [[ -e /mnt/disks/cromwell_root/gcs_delocalization.sh ]]
    	then
     		 CROMWELL_ROOT=/mnt/disks/cromwell_root
     	else
      		echo "Could not find Cromwell root under /cromwell_root (PAPI v2) or /mnt/disks/cromwell_root (GCP Batch), exiting."
      		exit 1
    	fi

    	echo $CROMWELL_ROOT
    }
    output {
        String cromwell_root_path = read_string(stdout())
    }
}