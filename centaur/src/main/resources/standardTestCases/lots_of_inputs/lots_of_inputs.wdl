task do_nothing {
  Array[File] f

  command {
    echo "no-op"
  }
  output {
    String o = read_string(stdout())
  }

  # We are localizing 10,000 tiny files, so the default disk size of 10 GB SSD [0] is NOT
  # going to cut it, having a measly IO per second rating of just 10 x 30 = 300 IOPS [1]
  # [0] https://cromwell.readthedocs.io/en/stable/RuntimeAttributes/
  # [1] https://cloud.google.com/compute/docs/disks/performance#zonal
  runtime {
    docker: "python:latest"
    disks: "local-disk 1000 SSD"
  }
}

task make_array {
  Int n
  command {
    python <<CODE
    for i in range(${n}):
      filename = 'file-' + str(i)
      with open(filename, 'w') as fp:
        fp.write(filename)
      print(filename)
    CODE
  }
  output {
    Array[File] a = glob("file-*")
  }
  runtime {
    docker: "python:latest"
    disks: "local-disk 1000 SSD"
  }
}

workflow lots_of_inputs {
  Int how_many_is_lots
  call make_array { input: n = how_many_is_lots }
  call do_nothing { input: f = make_array.a }

  output {
    Int out_count = length(make_array.a)
    String nothing_out = do_nothing.o
  }
}
