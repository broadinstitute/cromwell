version 1.0

workflow gpu_cuda_image {

    input {
      Array[String] driver_versions = [ "418.87.00" ]
    }

    scatter (driver_version in driver_versions) {
        call get_machine_info { input: driver_version = driver_version }
    }

    output {
      Array[String] modprobe_check = get_machine_info.modprobe_check
      Array[String] smi_check = get_machine_info.smi_check

      Array[File] modprobe_contents = get_machine_info.modprobe_content
      Array[File] smi_contents = get_machine_info.smi_content
    }
}

task get_machine_info {
    input {
        String driver_version
    }

	command <<<
		nvidia-modprobe --version > modprobe
		cat modprobe | grep -q "~{driver_version}" && echo "good" > modprobe_check || echo "bad" > modprobe_check
		nvidia-smi > smi
		cat smi | grep -q "~{driver_version}" && echo "good" > smi_check || echo "bad" > smi_check
	 >>>

    runtime {
        docker: "nvidia/cuda:9.0-cudnn7-devel-ubuntu16.04"
        bootDiskSizeGb: 20
        gpuType: "nvidia-tesla-k80"
        gpuCount: 1
        nvidiaDriverVersion: driver_version
        zones: "us-central1-c"
    }

    output {
      String modprobe_check = read_string("modprobe_check")
      String smi_check = read_string("smi_check")
      File modprobe_content = "modprobe"
      File smi_content = "smi"
    }
}
