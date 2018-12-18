version development

struct Person {
    String name
    Int age
}

workflow call_assignment_scattered_conditionals {
    Person harry = object {name: "Harry", age: 11}
    Array[Int] delays = range(5)

#    scatter(delay in delays) {
        Int delay = delays[1]
        if (delay % 2 == 1) {
            Int i = delay + 5
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
#    }

    output {
#        Array[Person?] harry_p = myTask
#        Array[Int?] is = i

        Person? herry_p = myTask
        Int? is = i

        Person? herry_p2 = myTask2

#        Array[Person] harry_p = myTask
#        Array[Int] is = i
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
