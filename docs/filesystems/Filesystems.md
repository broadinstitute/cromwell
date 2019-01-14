# Filesystems

Most workflows represent their inputs and outputs in the form of files. Those files are stored in filesystems. There exists many filesystems. This section describes which filesystems Cromwell supports.

## Overview

Filesystems are configurable. The `reference.conf`, which is the configuration inherited by any Cromwell instance, contains the following:

```hocon
# Filesystems available in this Crowmell instance
# They can be enabled individually in the engine.filesystems stanza and in the config.filesystems stanza of backends
# There is a default built-in local filesytem that can also be referenced as "local" as well.
filesystems {
  drs {
      class = "cromwell.filesystems.drs.DrsPathBuilderFactory"
      # Use to share a unique global object across all instances of the factory
      global {
        # Class to instantiate and propagate to all factories. Takes a single typesafe config argument
        class = "cromwell.filesystems.drs.DrsFileSystemConfig"
        config {
          martha {
            url = "https://martha-url-here"
            request.json-template = """{"url": "${drsPath}"}"""
          }
        }
      }
   }
  gcs {
    class = "cromwell.filesystems.gcs.GcsPathBuilderFactory"
  }
  oss {
    class = "cromwell.filesystems.oss.OssPathBuilderFactory"
  }
  s3 {
    class = "cromwell.filesystems.s3.S3PathBuilderFactory"
  }
  http {
    class = "cromwell.filesystems.http.HttpPathBuilderFactory"
  }
}
```

It defines the filesystems that can be accessed by Cromwell.
Those filesystems can be referenced by their name (`drs`, `gcs`, `oss`, `s3`, `http` and `local`) in other parts of the configuration.

**Note:**
- **OSS and S3 filesystems are experimental.** 
- **DRS filesystem has initial support only. Also, currently it works only with [GCS filesystem](../GoogleCloudStorage) in [PapiV2 backend](http://cromwell.readthedocs.io/en/develop/backends/Google).**


Also note that the local filesystem (the one on which Cromwell runs on) is implicitly accessible but can be disabled. 
To do so, add the following to any `filesystems` stanza in which the local filesystem should be disabled: `local.enabled: false`.

### Engine Filesystems

Cromwell is conceptually divided in an engine part and a backend part. One Cromwell instance corresponds to an "engine" but can have multiple backends configured.
The `engine.filesystems` section configures filesystems that Cromwell can use when it needs to interact with files outside of the context of a backend.

For instance, consider the following WDL:

```wdl
version 1.0

workflow my_workflow {
    String s = read_string("/Users/me/my_file.txt")
    output {
        String out = s
    }
}
```

This workflow is valid WDL and does not involve any backend, or even a task. However it does involve interacting with a filesystem to retrieve the content of `my_file.txt`
With a default configuration Cromwell will be able to run this workflow because the local filesystem is enabled by default.
If the file is located on a different filesystem (a cloud filesystem for instance), we would need to modify the configuration to tell Cromwell how to interact with this filesystem:

```hocon
engine {
  filesystems {
    gcs {
      auth = "application-default"
    }
  }
}
```

(See the [Google section](../backends/Google.md) for information about the `auth` field.)

We can now run this workflow

```wdl
version 1.0

workflow my_workflow {
    String s = read_string("gs://mybucket/my_file.txt")
    output {
        String out = s
    }
}
```

#### Default "engine" Filesystems

If you don't change anything in your own configuration file, the following default is inherited from `reference.conf`:
```
engine {
  filesystems {
    local {
      enabled: true
    }
    http {
      enabled: true
    }
  }
}
```

**Note**: since our configuration files are HOCON, to disable filesystems you *must* add `enabled: false` into your 
overriding configuration file. It is **not** sufficient to simply omit a filesystem from your stanza. 

For example: adding this to your configuration file will remove the `http` filesystem and leave `local` for use in the 
engine:
```
engine {
  filesystems {
    http {
      enabled: false
    }
  }
}
```

Whereas this example will leave `http` unchanged and merely re-assert the default enabling of `local`. In other
words, **this will do nothing**:
```
engine {
  filesystems {
    local {
      enabled: true
    }
  }
}
```

#### Engine filesystems and CWL

Note that CWL *always* needs to access file attributes from within the engine, so if you are using CWL, please make sure
that **every filesystem you might use** is added to the `engine.filesystems` stanza.

### Backend Filesystems

Similarly to the engine, you can also configure backend filesystems individually. Some backends might require the use of a specific filesystem.
For example, the [Pipelines API](../tutorials/PipelinesApi101.md) backend requires Google Cloud Storage.
Let's take another example:

```wdl
version 1.0

task my_pipelines_task {
    input {
        File input_file
    }
    String content = read_string(input_file)
    
    command {
        echo ~{content}
    }
    
    runtime {
        docker: "ubuntu"
    }
}
workflow my_workflow {
    call my_pipelines_task { input: input_file = "gs://mybucket/my_file.txt" }
}
```

Suppose this workflow is submitted to a Cromwell running a Pipelines API backend. This time the `read_string` function is in the context of a task run by the backend.
The filesystem configuration used will be the one in the `config` section of the Pipelines API backend.

### Supported Filesystems

-  Shared File System (SFS)

-  Google Cloud Storage (GCS) - [Cromwell Doc](GoogleCloudStorage.md) / [Google Doc](https://cloud.google.com/storage/)

-  Simple Storage Service (S3) - [Amazon Doc](https://aws.amazon.com/documentation/s3/)

-  Object Storage Service (OSS) - [Alibaba Cloud Doc](https://www.alibabacloud.com/product/oss)

-  HTTP - support for `http` or `https` URLs for [workflow inputs only](http://cromwell.readthedocs.io/en/develop/filesystems/HTTP)

- File Transfer Protocol (FTP) - [Cromwell Doc](FileTransferProtocol.md)
