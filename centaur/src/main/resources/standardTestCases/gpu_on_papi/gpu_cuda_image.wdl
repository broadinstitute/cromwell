version 1.0

workflow gpu_cuda_image {

    call get_machine_info

    output {
      String smi_check = get_machine_info.smi_check
      File smi_contents = get_machine_info.smi_content
    }
}

task get_machine_info {

	command <<<
		nvidia-smi > smi
		cat smi | grep -q "Tesla T4" && echo "gpu_good" > smi_check || echo "bad" > smi_check
		cat smi | grep -q "15360MiB" && echo "vram_good" >> smi_check || echo "bad" >> smi_check
	 >>>

    runtime {
        docker: "nvidia/cuda:12.6.2-cudnn-devel-ubuntu24.04"
        bootDiskSizeGb: 20
        gpuType: "nvidia-tesla-t4"
        gpuCount: 1
        zones: "us-central1-c"
    }

    output {
      String smi_check = read_string("smi_check")
      File smi_content = "smi"
    }
}
