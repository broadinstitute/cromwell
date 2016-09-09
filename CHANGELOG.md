# Cromwell Change Log

## 0.21

* Add support for Google Private IPs through `noAddress` runtime attribute. If set to true, the VM will NOT be provided with a public IP address.
*Important*: Your project must be whitelisted in "Google Access for Private IPs Early Access Program". If it's not whitelisted and you set this attribute to true, the task will hang.
  Defaults to `false`.
  e.g:
  
```
task {
    command {
        echo "I'm private !"
    }
    
    runtime {
        docker: "ubuntu:latest"
        noAddress: true
    }
}
```