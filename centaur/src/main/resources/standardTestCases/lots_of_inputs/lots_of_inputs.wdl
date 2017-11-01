task do_nothing {
  Array[File] f

  command {
    echo "no-op"
  }
  output {
    String o = read_string(stdout())
  }
  runtime {
    docker: "python:latest"
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
  }
}

workflow lots_of_inputs {
  call make_array {input: n=400}
  call do_nothing {input: f=make_array.a}
}
