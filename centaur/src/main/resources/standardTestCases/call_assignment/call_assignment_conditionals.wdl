version development

struct Person {
    String name
    Int age
}

workflow call_assignment_conditionals {
    Person harry = object {name: "Harry", age: 11}
    Int delay = 1
    if (delay % 2 == 1) {
        call myTask { input:
            a = harry,
            delay = delay
        }
    }
    if (delay % 2 == 0) {
        call myTask as myTask2 { input:
            a = harry,
            delay = delay
        }
    }

    output {
        Person? herry_p = myTask
        Person? herry_p2 = myTask2
    }
}

task myTask {
    input {
        Person a
        Int delay
    }
    command <<<
        sleep ~{delay * 2}
        echo "hello my name is ~{a.name} and I am ~{a.age} years old"
    >>>
    output {
        String name = a.name + " Potter ~{delay}"
        Int age = a.age + delay
    }
}
