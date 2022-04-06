MINIO_ACCESS_KEY=minioadmin MINIO_SECRET_KEY=minioadmin minio server test-data 



curl -L -O https://github.com/ohsu-comp-bio/funnel/releases/download/0.10.1/funnel-linux-amd64-0.10.1.tar.gz

tar xvzf funnel-linux-amd64-0.10.1.tar.gz funnel



GenericS3:
  - Disabled: false
    Endpoint: "192.168.86.37:9000"
    Key: "minioadmin"
    Secret: "minioadmin"



```
{
  "name": "Hello world",
  "inputs": [{
    # URL to download file from.
    "url": "s3://test/README.md",
    # Path to download file to.
    "path": "/inputs/hello.txt"
  }],
  "outputs": [{
    # URL to upload file to.
    "url": "s3://test/output.txt",
    # Local path to upload file from.
    "path": "/outputs/stdout"
  }],
  "executors": [{
      # Container image name.
      "image": "alpine",
      # Command to run (argv).
      "command": ["cat", "/inputs/hello.txt"],
      # Capture the stdout of the command to /outputs/stdout
      "stdout": "/outputs/stdout"
  }]
}
```

./funnel task create funnel-test.json 



sbt assembly

```
minio {
  application-name = "cromwell"
  auths = [
    {
      name = "default"
      endpoint = "http://localhost:9000"
      access-key = "minioadmin"
      secret-key = "minioadmin"
    }
  ]
}

engine {
  filesystems {
    minio {
      auth = "default"
    }
  }
}

backend {
  providers {
    Local {
        actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"
        config {
          concurrent-job-limit = 12
        }
    }
  }
}

filesystems { 
  minio {
      class = "cromwell.filesystems.minio.MinioPathBuilderFactory"
  }
}
```

```
# include required(classpath("application"))

minio {
  application-name = "cromwell"
  auths = [
    {
      name = "default"
      endpoint = "http://localhost:9000"
      access-key = "minioadmin"
      secret-key = "minioadmin"
    }
  ]
}

engine {
  filesystems {
    minio {
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
        root = "cromwell-executions"
        dockerRoot = "/cromwell-executions"
        endpoint = "http://127.0.0.1:8000/v1/tasks"
        default-runtime-attributes {
          cpu: 1
          failOnStderr: false
          continueOnReturnCode: 0
          memory: "2 GB"
          disk: "2 GB"
          preemptible: false
        }
      }
    }
  }
}

filesystems { 
  minio {
      class = "cromwell.filesystems.minio.MinioPathBuilderFactory"
  }
}
```


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

```
{
    "testwf.wc.in_file" : "minio://test/README.md"
}
```

java -Dconfig.file=cromwell_funnel.conf -jar ./server/target/scala-2.12/cromwell-79-2802c31-SNAP.jar run test.wdl -i test.json 
