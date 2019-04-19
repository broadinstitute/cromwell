
task A {
  command {
    echo "Enfin un peu de francais pour contrer ce raz-de-marée anglais !" > out
    echo "Jacques Chirac fait du jetski sur la Seine en costume traditionnel russe" > out2
  }
  output {
    File out = "out"
    File out2 = "out2"
  }
  runtime {
    docker: "ubuntu:latest"
  }
}

task B {
  command {
     echo "Je contre avec un bonnet peruvien et tire une carte chance" > out
     echo "Kamoulox !" > out2
  }
  output {
     Array[File] outs = ["out", "out2"]
  }
  runtime {
    docker: "ubuntu:latest"
  }
}

task C {
  command {
    cat > out <<END
    (_)
     !_________________________________________
     !*  *  *  *  * |##########################|
     ! *  *  *  *  *|                          |
     !*  *  *  *  * |##########################|
     ! *  *  *  *  *|                          |
     !*  *  *  *  * |##########################|
     ! *  *  *  *  *|                          |
     !*  *  *  *  * |##########################|
     !~~~~~~~~~~~~~~~                          |
     !#########################################|
     !                                         |
     !#########################################|
     !                                         |
     !###################################JGS###|
     !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     !
     !
     !
     !
     !
     !
     !
    END
  }
  output {
      File out = "out"
  }
  runtime {
    docker: "ubuntu:latest"
  }
}

workflow wfoutputs {
  call A
  call B
  call C
  output {
    A.*
    B.outs
  }
}
