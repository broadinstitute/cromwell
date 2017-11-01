task make_file {
  command {
    echo There once was a fisher named Fisher > out
    echo who fished for some fish in a fissure >> out
    echo Till a fish with a grin >> out
    echo pulled the fisherman in >> out
    echo Now theyre fishing the fissure for Fisher >> out
  }
  runtime {
    docker: "ubuntu:latest"
  }
  output {
    File out = "out"
  }
}

task maybe_cats {
  Pair[File?, File?] opt_file_pair

  command {
    touch left
    touch right
    ${ "cat " + opt_file_pair.left + " > left"}
    ${ "cat " + opt_file_pair.right + " > right"}
  }
  runtime {
    docker: "ubuntu:latest"
  }
  output {
    Pair[String, String] result = (read_string("left"), read_string("right"))
  }
}

task array_cat {
  Array[File] file_arr

  command {
    cat ${sep=" " file_arr}
  }
  runtime {
    docker: "ubuntu:latest"
  }
  output {
    String result = read_string(stdout())
  }
}

task map_cat {
  Map[String, File?] file_map

  command {
    ${ "cat " + file_map["always"] }
    ${ "cat " + file_map["truth"] }
    ${ "cat " + file_map["untruth"] }
  }
  runtime {
    docker: "ubuntu:latest"
  }
  output {
    String result = read_string(stdout())
  }
}

workflow complex_types_files {
  call make_file
  if (true) { call make_file as make_true }
  if (false) { call make_file as make_false }

  call maybe_cats { input: opt_file_pair = (make_true.out, make_false.out) }
  Array[File] file_arr = select_all([make_file.out, make_true.out, make_false.out])
  call array_cat { input: file_arr = file_arr }
  Map[String, File?] file_map = { "always": make_file.out, "truth": make_true.out, "untruth": make_false.out }
  call map_cat { input: file_map = file_map }

  output {
    Array[String] result = [ maybe_cats.result.left, maybe_cats.result.right, array_cat.result, map_cat.result ]
  }
}
