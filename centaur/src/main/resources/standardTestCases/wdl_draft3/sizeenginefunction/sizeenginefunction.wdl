version 1.0

task make_a_file {
  command {
    echo "hey look here's a file" > a_file
  }
  output {
    File a_file = "a_file"
  }
  runtime {docker: "ubuntu:latest"}
}

task size {
  input {
    File? unsupplied
    File input_file
  }

  command <<<
    echo "this file is 22 bytes" > created_file
    ~{"cat " + input_file}
    ~{ if defined(unsupplied) then "cat " + unsupplied else "" }
  >>>

  output {
    Float input_file_size = size(input_file)
    Float input_file_size2 = size(stdout())
    Float unsupplied_size = size(unsupplied)
    Int total_input_size = round(size(unsupplied) + size(input_file))
    Int created_file_size = round(size("created_file"))
    Float created_file_size_in_k = size("created_file", "K")
    Float created_file_size_in_kb = size("created_file", "KB")
    Float created_file_size_in_m = size("created_file", "M")
    Float created_file_size_in_mb = size("created_file", "MB")
    Float created_file_size_in_g = size("created_file", "G")
    Float created_file_size_in_gb = size("created_file", "GB")
    Float created_file_size_in_t = size("created_file", "T")
    Float created_file_size_in_tb = size("created_file", "TB")
    Float created_file_size_in_ki = size("created_file", "Ki")
    Float created_file_size_in_kib = size("created_file", "KiB")
    Float created_file_size_in_mi = size("created_file", "Mi")
    Float created_file_size_in_mib = size("created_file", "MiB")
    Float created_file_size_in_gi = size("created_file", "Gi")
    Float created_file_size_in_gib = size("created_file", "GiB")
    Float created_file_size_in_ti = size("created_file", "Ti")
    Float created_file_size_in_tib = size("created_file", "TiB")
  }

  runtime {docker: "ubuntu:latest"}
}

workflow sizeenginefunction {
  scatter(x in range(3)) {
    call make_a_file
  }

  call size {
      input: input_file = make_a_file.a_file[0]
  }

  scatter(made in make_a_file.a_file) {
    Int workflow_size = round(size(made))
  }

  output {
    Array[Int] workflow_sizes = workflow_size
    Float input_file_size = size.input_file_size
    Float input_file_size2 = size.input_file_size2
    Float unsupplied_size = size.unsupplied_size
    Int total_input_size = size.total_input_size
    Int created_file_size = size.created_file_size
    Float created_file_size_in_k = size.created_file_size_in_k
    Float created_file_size_in_kb = size.created_file_size_in_kb
    Float created_file_size_in_m = size.created_file_size_in_m
    Float created_file_size_in_mb = size.created_file_size_in_mb
    Float created_file_size_in_g = size.created_file_size_in_g
    Float created_file_size_in_gb = size.created_file_size_in_gb
    Float created_file_size_in_t = size.created_file_size_in_t
    Float created_file_size_in_tb = size.created_file_size_in_tb
    Float created_file_size_in_ki = size.created_file_size_in_ki
    Float created_file_size_in_kib = size.created_file_size_in_kib
    Float created_file_size_in_mi = size.created_file_size_in_mi
    Float created_file_size_in_mib = size.created_file_size_in_mib
    Float created_file_size_in_gi = size.created_file_size_in_gi
    Float created_file_size_in_gib = size.created_file_size_in_gib
    Float created_file_size_in_ti = size.created_file_size_in_ti
    Float created_file_size_in_tib = size.created_file_size_in_tib
  }
}
