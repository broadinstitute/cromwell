task t {
    String i
    command {
        echo "${i}"
    }
    runtime {
        docker: "ubuntu:latest"
    }
    output {
        String o = read_string(stdout())
    }
}

workflow declarations_as_nodes {
    call t as t1 { input: i = "hello" }         # hello
    
    String a = t1.o + " world"                  # hello world
    
    call t as t2 { input: i = a }               # hello world
    
    Array[String] arr = [t1.o, t2.o]            # [hello, hello world]
    
    scatter(i in arr) {
        call t as t3 { input: i = i }           # hello | hello world
        String b = i + t3.o                     # hellohello | hello worldhello world
        call t as t4 { input: i = b }           # hellohello | hello worldhello world
        String c = t3.o + " " + t4.o            # hello hellohello | hello world hello worldhello world
    }
    
    Array[String] d = c                         # [hello hellohello, hello world hello worldhello world]
     
    output {
        String o1 = a                           # hello world
        Array[String] o2 = b                    # [hellohello, hello worldhello world]
        Array[String] o3 = c                    # [hello hellohello, hello world hello worldhello world]
        Array[String] o4 = d                    # [hello hellohello, hello world hello worldhello world]
    }
}
