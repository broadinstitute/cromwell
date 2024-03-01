version development-1.1

struct Animal {
    String name
    Boolean isGood
}

task create_dog {
    input {
        String name_input
    }

    command { echo "all dogs are good" }
    output {
        Animal dog = Animal{name: name_input, isGood: true}
    }
}

workflow struct_literal {
    call create_dog { input: name_input = "doggo" }
    output {
        String name_output = create_dog.dog.name
    }
}
