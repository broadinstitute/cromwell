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
        Plant standard_plant_input = Plant{color: "green", tasty: true}
        Animal standard_animal_input = Animal{name: "mittens", isGood: false}
        Animal omitted_optional_animal = Animal{name: "boots"}
    }

    meta {
        volatile: true
    }

    runtime {
        docker: "ubuntu:latest"
    }

    command { echo "all dogs are good" }

    output {
        Plant standard_plant_forwarded = standard_plant_input
        Animal standard_animal_forwarded = standard_animal_input
        Plant plant_output_literal = Plant{color: "red", tasty: true}
        Animal animal_output_literal = Animal{name: "Drainus", isGood: true}
        Animal omitted_output_forwarded_animal = omitted_optional_animal
        Animal omitted_output_literal_animal = Animal{name: "Alfonso"}
    }
}

workflow struct_literal {
    call test_struct_parsing
    output {
        Plant forwarded_input_1 = test_struct_parsing.standard_plant_forwarded
        Animal forwarded_input_2 = test_struct_parsing.standard_animal_forwarded
        Plant literal_output_1 = test_struct_parsing.plant_output_literal
        Animal literal_output_2  = test_struct_parsing.animal_output_literal
        Animal omitted_optional_1 = test_struct_parsing.omitted_output_forwarded_animal
        Animal omitted_optional_2 = test_struct_parsing.omitted_output_literal_animal
    }
}
