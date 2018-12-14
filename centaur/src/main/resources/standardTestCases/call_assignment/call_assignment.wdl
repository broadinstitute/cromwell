version development

struct Person {
    String name
    Int age
}

workflow call_assignment {
    Person harry = object {name: "Harry", age: 11}
    call myTask {
        input:
            a = harry
    }
    output {
      Person harry_p = myTask
    }
}

task myTask {
    input {
      Person a
    }
    command <<<
        echo "hello my name is ~{a.name} and I am ~{a.age} years old"
    >>>
    output {
      String name = a.name + " Potter"
      Int age = a.age * 2
    }
}
