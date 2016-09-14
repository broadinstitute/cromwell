# Cromwell Change Log

## 0.20

* The default per-upload bytes size for GCS is now the minumum 256K
instead of 64M. There is also an undocumented config key
`google.upload-buffer-bytes` that allows adjusting this internal value.

* Updated Docker Hub hash retriever to parse json with [custom media
types](https://github.com/docker/distribution/blob/05b0ab0/docs/spec/manifest-v2-1.md).

* Added a `/batch` submit endpoint that accepts a single wdl with
multiple input files.

* The `/query` endpoint now supports querying by `id`, and submitting
parameters as a HTTP POST.

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

* The Local and the SGE backend have been merged into a generic
Shared File System (SFS) backend. This updated backend can be configured
to work with various other command line dispatchers such as LSF. See the
[README](README.md#sun-gridengine-backend) for more info.

* On the JES and SFS backends, task `command` blocks are now always
passed absolute paths for input `File`s.

* On the SFS backends, the call directory now contains two sub-directories:
    * `inputs` contains all the input files that have been localized for this task (see next below for more details)
    * `execution` contains all other files (script, logs, rc, potential outputs etc...)
    
* On the SFS backends, inputs are now localized as follow:
    The absolute path of the directory containing the input is hashed.
    This hash is used to create a directory under `inputs` where the file will be localized
    e.g:
    
    ```
    ├── execution
    │   ├── rc
    │   ├── script
    │   ├── script.background
    │   ├── script.submit
    │   ├── stderr
    │   ├── stderr.background
    │   ├── stdout
    │   └── stdout.background
    └── inputs
        ├── 2c142eed1d519eabc0520becff443327
        │   ├── file3.txt
        │   └── file4.txt
        └── 63a9f0ea7bb98050796b649e85481845
            ├── file1.txt
            └── file2.txt
    ```
    
    where the file*.txt input files path structure would be similar too:
    
    ```
    root
    ├── file1.txt       -> will be localized under hash("root")/file1.txt
    ├── file2.txt       -> will be localized under hash("root")/file2.txt
    └── mydir
        ├── file3.txt   -> will be localized under hash("root/mydir")/file3.txt
        └── file4.txt   -> will be localized under hash("root/mydir")/file4.txt
    ```