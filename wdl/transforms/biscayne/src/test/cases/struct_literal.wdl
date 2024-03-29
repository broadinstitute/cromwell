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
        Plant standard_plant_input
        Animal standard_animal_input
    }

    runtime {
        docker: "ubuntu:latest"
    }

    command {
       echo "all dogs are good"
    }

    output {
        Plant standard_plant_forwarded = standard_plant_input
        Animal standard_animal_forwarded = standard_animal_input
        Plant plant_output_literal = Plant{color: "red", tasty: true}
    }
}

workflow struct_literal {
    call test_struct_parsing {
         input: standard_plant_input = Plant{color: "green", tasty: true},
            standard_animal_input = Animal{name: "mittens", isGood: false}
         }
}
