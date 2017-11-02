task sleep_exit {
    Int s
    Boolean flow_control = true
    Int exit_code = 0
    command {
        sleep ${s}
        exit ${exit_code}
    }
    runtime {
       docker: "ubuntu:latest"
    }
    output {
        Boolean done = true
    }
}

workflow restart_while_failing {
    # A will fail after 2 minutes
    call sleep_exit as A { input: s = 120, exit_code = 1 }
    
    # B starts at the same time as A and finishes immediately
    call sleep_exit as B { input: s = 10 }
    # B1 starts as soon as B is done
    call sleep_exit as B1 { input: s = 180, flow_control = B.done }
    # B2 depends on B1. Note that B2 should never be run because by the time we get to it, A will have failed.
    call sleep_exit as B2 { input: s = 10, flow_control = B1.done }
    
    Array[Int] sleep_array = [1, 60, 20, 60, 30]
    # Independently, we start a scatter with a few jobs, some of them will still be running when we restart 
    scatter(i in sleep_array) {
       call sleep_exit { input: s = i }
    }
    
    # Cromwell will restart as soon as B1 is running
    # The goal is to make sure that all jobs in this workflow reach a terminal status after the Cromwell restart,
    # and the workflow itself is failed
}
