version development-1.1


struct Animal {
 String name
 Boolean? isGood
}

task create_dog {
    input {
        String name_input
    }
    runtime {
        docker: "ubuntu:latest"
    }

    command { echo "all dogs are good" }
    output {
        Animal cat = Animal{name: "mittens", isGood: false}
        Animal dog = Animal{name: name_input, isGood: true}
    }
}

workflow struct_literal {
    call create_dog { input: name_input = "doggo" }
    output {
        Boolean? areDogsGood = create_dog.dog.isGood
        Boolean? areCatsGood = create_dog.cat.isGood
        String name_output = create_dog.dog.name
    }
}
