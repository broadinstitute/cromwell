version development-1.1

struct Plant {
    String color
    Boolean tasty
}

struct Animal {
    String name
    Boolean? isGood
}

task test_struct_parsing {
    input {
        Plant p1 = Plant{color: "green", tasty: true}
        Plant p2 = Plant{color: "brown", tasty: false}
        Plant p3 = Plant{color: "red", tasty: false}
        Animal a1 = Animal{name: "mittens", isGood: false}
    }

    meta {
        volatile: true
    }

    runtime {
        docker: "ubuntu:latest"
    }

    command { echo "all dogs are good" }

    output {
        Plant o1 = Plant{color: "green", tasty: true}
        Plant o2 = Plant{color: "green", tasty: true}
        Plant o3 = Plant{color: "green", tasty: true}
        Animal o4 = Animal{name: "BlobStorage", isGood: true}
    }
}

workflow struct_literal {
    call test_struct_parsing
    output {
        Boolean tasty = test_struct_parsing.o1.tasty
        String pet = test_struct_parsing.o4.name
    }
}
