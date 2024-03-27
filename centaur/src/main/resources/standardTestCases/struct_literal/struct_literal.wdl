version development-1.1

struct Plant {
  String color
  Int id
}

struct Fungi {
    File fungiFile
}


struct Animal {
  Plant jacket
  Fungi hat
}

task a {
    input {
      Plant in_plant_literal = Plant{color: "red", id: 44}
    }

  command {
    echo "${in_plant_literal.id}"
  }

  output {
    Animal out_animal = Animal{jacket: Plant{color: "green", id: 10}, hat: Fungi{fungiFile: stdout()}}
  }

  runtime {
    docker: "ubuntu:latest"
 }

  meta {
    volatile: true
  }
}

task b {
  input {
  Animal in_animal
  }

  command {
    cat ${in_animal.hat.fungiFile}
  }

  output {
    Int out = read_int(stdout())
  }

  runtime {
   docker: "ubuntu:latest"
 }

  meta {
   volatile: true
  }
}

workflow struct_literal {
  call a
  call b {input: in_animal=a.out_animal}
  output {
    Int out = b.out
  }
}
