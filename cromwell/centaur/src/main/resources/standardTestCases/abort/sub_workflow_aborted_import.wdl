task aborted {
    Boolean go
    # This sleep is long on purpose: we'll check that the VM doesn't exist anymore as a proof that the job
    # was successfully aborted, but we don't want it to die on its own before that because the job ended
    command {
        echo "Abort incoming"
        sleep 1200
    }
    runtime {
       docker: "ubuntu:latest"
    }
    output {
      Boolean done = true
   }
}
workflow inner_abort {
  Boolean go
  
  call aborted { input: go = go }
  
  output {
    Boolean done = aborted.done
  }
}
