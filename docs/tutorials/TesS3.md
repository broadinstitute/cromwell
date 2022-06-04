
## Getting started with TES and Generic S3

The TES API allows Cromwell to plug into a number of different job scheduling systems.
Combined with object storage, this allows jobs to be easily moved to new computing environments.
S3, while primarily utilized by AWS, can also be provided by backends, including CEPH and MinIO.
Because these systems provide the S3 protocol, but don't use other AWS systems, the driver is called
`GenericS3`. For this tutorial, we will setup a private S3 instance, connect it to a TES server and 
configure Cromwell to utilize these systems.

### Setting up a custom S3 instance on a private server.

(MinIO)[https://min.io] is a lightweight S3 server. To download and run a copy, visit the instructions 
at https://min.io/download.

For linux:
```
wget https://dl.min.io/server/minio/release/linux-amd64/minio
chmod +x minio
```

To launch a test instance with the default access and secret keys. This will turn the directory `data-dir` into 
an S3 storage point:
```
MINIO_ROOT_USER=minioadmin MINIO_ROOT_PASSWORD=minioadmin ./minio server data-dir --console-address ":9001"
```

This will print several messages, one will look like:
```
API: http://192.168.86.37:9000  http://172.22.0.1:9000  http://172.18.0.1:9000  http://172.19.0.1:9000  http://172.20.0.1:9000  http://172.21.0.1:9000  http://172.17.0.1:9000  http://127.0.0.1:9000                     
RootUser: minioadmin 
RootPass: minioadmin 
```

Please note the API address URL, this will be used later on to configure file access in Funnel and Cromwell.
You should also see a URL for the web console. Visit that website and create a bucket named `cromwell`.
Alternatively, because Minio reflects the contents of the `test-data` directory, simply run the command:

```
mkdir test-data/cromwell
```

Also upload a file to use as an input to workflows, called `input.data`. Again, you can also just copy 
a file into the `test-data` directory:
```
cp README.md test-data/cromwell/input.data
```

This file will now have the address `s3://cromwell/input.data`


### Start local TES server

Download a copy of Funnel:
```
curl -L -O https://github.com/ohsu-comp-bio/funnel/releases/download/0.10.1/funnel-linux-amd64-0.10.1.tar.gz
tar xvzf funnel-linux-amd64-0.10.1.tar.gz funnel
```

Create a configuration file `funnel.config` that directs it to use the S3 storage server you just set up.
This will only run jobs in local mode, on the current machine. For more advanced configurations that use engines like
SLURM or K8S, see the Funnel manual.
```
GenericS3:
  - Disabled: false
    Endpoint: "192.168.86.37:9000"
    Key: "minioadmin"
    Secret: "minioadmin"
```

Start the Funnel server:
```
./funnel server run -c funnel.config 
```

You can run a test job, to make sure that Funnel can access the files and execute jobs.
Create the file `funnel-test.json`:
```
{
  "name": "Hello world",
  "inputs": [{
    "url": "s3://cromwell/input.data",
    "path": "/inputs/hello.txt"
  }],
  "outputs": [{
    "url": "s3://cromwell/output.txt",
    "path": "/outputs/stdout"
  }],
  "executors": [{
      "image": "alpine",
      "command": ["cat", "/inputs/hello.txt"],
      "stdout": "/outputs/stdout"
  }]
}
```
Then submit it:
```
./funnel task create funnel-test.json 
```
The file `s3://cromwell/output.txt` should now exist. 


### Connecting Cromwell

The configuration connect Cromwell to the MinIO server and Funnel:

```
# include required(classpath("application"))

genericS3 {
  name = "default"
  endpoint = "http://192.168.86.37:9000"
  access-key = "minioadmin"
  secret-key = "minioadmin"
}

filesystems { 
  genericS3 {
      class = "cromwell.filesystems.s3.GenericS3PathBuilderFactory"
  }
}

engine {
  filesystems {
    genericS3 {
      auth = "default"
    }
  }
}

backend {
  default = TES

  providers {
    TES {
      actor-factory = "cromwell.backend.impl.tes.TesBackendLifecycleActorFactory"
      config {
        root = "s3://cromwell/executions/"
        dockerRoot = "/docker-cromwell-executions"
        endpoint = "http://192.168.86.37:8000/v1/tasks" # replace with your TES API endpoint
        default-runtime-attributes {
          cpu: 1
          failOnStderr: false
          continueOnReturnCode: 0
          memory: "2 GB"
          disk: "2 GB"
          preemptible: false
        }

        filesystems {
          genericS3 {
            auth = "default"
          }
        }
      }
    }
  }
}
```

### Running a workflow

Using the following workflow `test.wdl`:
```
workflow testwf {
  call wc
}

task wc {
  File in_file
  command {
    cat ${in_file} | wc -l
  }
  output {
    Int count = read_int(stdout())
  }
}
```

The the inputfile `test.json`:
```
{
    "testwf.wc.in_file" : "s3://cromwell/input.data"
}
```

Run the command:
```
java -Dconfig.file=cromwell_funnel.conf -jar cromwell.jar run test.wdl -i test.json 
```

You should now see the contents of the workflow under `s3://cromwell/executions`
