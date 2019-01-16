task A {
  command {
    echo "The creatures outside looked from pig to man, and from man to pig, and from pig to man again: but already it was impossible to say which was which." > B1
    echo "But it was all right, everything was all right, the struggle was finished. He had won the victory over himself. He loved Big Brother." > B2
  }
  output {
    Array[File] outs = [ "B1", "B2" ]
  }
  runtime {
    docker: "ubuntu:latest"
  }
}

workflow wfoutputs {
  Array[Int] xs = [1,2,3,4,5,6,7,8,9,10]

  scatter ( x in xs ) {
    call A
  }
}
