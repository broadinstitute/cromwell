**Alibaba Cloud BCS Backend**

This backend adds support for execution jobs on Alibaba Cloud's BatchCompute service in a workflow.

### Configuring Backend

The backend is specified via the actor factory `BcsBackendLifecycleActorFactory`:

```hocon
backend {
  providers {
    BCS {
      config {
        actor-factory = "cromwell.backend.impl.bcs.BcsBackendLifecycleActorFactory"
        # ... other configuration
      }
    }
  }
}
```

You'll likely also want to change the default backend to your new backend, by setting this configuration value:

```hocon
backend {
  providers {
    default = BCS
  }
}
```

Before reading further in this section please see the [Getting started on Alibaba Cloud](../tutorials/BCSIntro.md) for instructions on configuring to Alibaba Cloud services.

The configuration file for Alibaba Cloud will look like the following.

```hocon
backend {
  providers {
    BCS {
      config {
        actor-factory = "cromwell.backend.impl.bcs.BcsBackendLifecycleActorFactory"
        root = "oss://<test-bucket>/cromwell-dir"
        region = "<test-region>"
        access-id = "<test-access-id>"
        access-key = "<test-access-key>"
       
        filesystems {
        # ... to be filled in
        }
        
        default-runtime-attributes {
        # ... to be filled in
        }
      }
    }
  }
}
```

- `<test-bucket>` : OSS bucket name.
- `<test-region>` : Region in Alibaba Cloud chosen to deploy cromwell, it must be the same as the region of `<test-bucket>`.
- `<test-access-id>` : Access ID to access Alibaba Cloud services through restful API.
- `<test-access-key>` : Access key to access Alibaba Cloud services through restful API.

The values above are necessary for Cromwell to submit and poll status of workflow jobs to and from Alibaba Cloud BatchCompute service.
The `filesystems` stanza in the backend config defines how to configure a filesystem in Alibaba Cloud. Details of filesystem related configurations will be explained in the next section. 

### File Systems

Currently, this backend only works with objects on an Alibaba Cloud OSS filesystem. It's necessary to supply all values in 
the configuration key `backend.providers.BCS.config.filesystems.auth` in order to read/write OSS file system objects in Alibaba Backend jobs. A typical config looks like this:

- `<test-oss-endpoint>` - API endpoint to access OSS bucket `<test-bucket>`.
- `<test-access-id>` - Access ID to access Alibaba Cloud services through restful API. 
- `<test-access-key>` - Access key to access Alibaba Cloud services through restful API. 
- `<refresh-interval>` - The interval of auth refreshing if you are using an STS(Alibaba Cloud Security Token Service) way to access the OSS filesystem.

```hocon
backend {
  providers {
    BCS {
      config {
        # BCS related configurations mentioned above
       
        filesystems {
          oss {
            auth {
              endpoint = "<test-oss-endpoint>"
              access-id = "<test-access-id>"
              access-key = "<test-access-key>"
            }
            refresh-interval = 1800
          }
        }
        
        default-runtime-attributes {
        # ... to be filled in
        }
      }
    }
  }
}
```

### Runtime Attributes

This backend supports additional runtime attributes that are specified in the configuration key `backend.providers.BCS.config.runtime-attributes`. 
It uses the same syntax as specifying runtime attributes in a task in WDL. A typical runtime attributes example for BCS backend looks like this:

```hocon
backend {
  providers {
    BCS {
      config {
        # BCS and OSS related configurations mentioned above
       
        default-runtime-attributes {
          cluster: "OnDemand ecs.sn1ne.large img-ubuntu-vpc"
          imageId: "img-ubuntu-vpc"
          mounts: "oss://<test-bucket>/inputs/ /home/inputs/ false"
          dockerTag: "ubuntu/latest oss://<test-bucket>/registry/ubuntu/"
          docker: "registry.cn-shanghai.aliyuncs.com/batchcompute/myubuntu:0.2"
          userData: "key value"
          reserveOnFail: true
          autoReleaseJob: true
          verbose: false
          systemDisk: "cloud 50"
          dataDisk: "cloud 250 /home/data/"
          timeout: 3000
        }
      }
    }
  }
}
```

#### cluster

There are two different ways of specifying an Alibaba Cloud BatchCompute cluster in which workflow jobs run.

- Reserved cluster - A pre-created cluster ID in BatchCompute service like this:

```hocon
      default-runtime-attributes {
        cluster: "cls-your-cluster-id"
      }
```

- Auto cluster - Cluster configuration to create a new runtime cluster bound to the workflow job:

  - `<resource-type>` - Type of resource, can only support `OnDemand` and `Spot` currently.
  - `<instance-type>` - Type of VM instance. Go to <a href="https://help.aliyun.com/document_detail/25378.html" target="_blank">Alibaba Cloud BatchCompute Instance Type</a> to choose a suitable type for you.
  - `<image-id>` - Image ID of Alibaba Cloud BatchCompute service to create a VM.

```hocon
      default-runtime-attributes {
        cluster: "<resource-type> <instance-type> <image-id>"
        # Maybe like cluster: "OnDemand ecs.sn1ne.large img-ubuntu"
      }
```

#### imageId

The BCS job image ID can be specified by the runtime `cluster`, while if you are going to use Call Caching, another optional runtime 
`imageId` should be specified. The change of image ID will lead to a cache miss and the call will be executed as normal.

```hocon
      default-runtime-attributes {
        imageId: "img-ubuntu-vpc"
      }
```

#### mounts

BCS jobs can mount both OSS and [Alibaba Cloud NAS](https://www.aliyun.com/product/nas) to local filesystem as a file or a directory in VM.
It uses distribute-caching and lazy-load techniques to optimize concurrently read requests of the OSS file system. 
You can mount your OSS objects to VM like this:

- `<mount-src>` - An OSS object path or OSS prefix or NAS address to mount from, such as
  `oss://<test-bucket>/inputs/ /home/inputs/ false` for OSS
  and `nas://0266e49fea-yio75.cn-beijing.nas.aliyuncs.com:/ /home/nas/ true` for NAS. See the [NAS mount](https://www.alibabacloud.com/help/doc-detail/50494.htm) for more details of NAS mount.
- `<mount-destination>` - An unix file path or directory path to mount to in VM.
- `<write-support>` - Writable for mount destination, only works for directory.

```hocon
default-runtime-attributes {
  mounts: "<mount-src> <mount-destination> <write-support>"
}
```



#### docker

This backend supports docker images pulled from OSS registry or Alibaba Cloud Container Registry.

##### OSS registry
```hocon
default-runtime-attributes {
  dockerTag: "<docker-image> <oss-registry-path>"
}
```

- `<docker-image>` - Docker image name such as: ubuntu:latest.
- `<oss-registry-path>` - Image path in OSS filesyetem where you pushed your docker image.


##### Alibaba Cloud Container Registry

```hocon
default-runtime-attributes {
  docker: "<docker-image-with-tag>"
}
```
- `docker-image-with-tag` - Docker image stored in Alibaba Cloud Container Registry, such as `registry.cn-shanghai.aliyuncs.com/batchcompute/myubuntu:0.2`.

#### userData

If a runtime cluster is specified, it's possible to pass some environment variables to VM when running BCS jobs.
It looks like this:

```hocon
 default-runtime-attributes {
   userData: "key1 value1, key2, value2"
 }
```

#### autoReleaseJob

The Alibaba Cloud BatchCompute service limits the number of simultaneous jobs per user. Jobs created by the backend are
deleted when the job finishes. However it is possible to tell the BCS backend to not delete the job when the call
finishes by setting `autoReleaseJob` to `false`:

```hocon
 default-runtime-attributes {
   autoReleaseJob: false
 }
```

#### systemDisk

If it's necessary to run a job with a particular system disk type or disk size, a runtime attribute named `systemDisk` can be used to
specify disk type and size.

- `<disk-type>` - Disk type to be used, can only support `cloud` or `cloud_efficiency` currently.
- `<disk-size-in-GB>` - Disk size to be used.

```hocon
 default-runtime-attributes {
   systemDisk: "<disk-type> <disk-size-in-GB>"
 }
```

#### dataDisk

The system disk size can support up to 500GB. One can mount another data disk in VM if needed.

- `<disk-type>` - Disk type to be used, can only support `cloud` or `cloud_efficiency` currently.
- `<disk-size-in-GB>` - Disk size to be used.
- `<mount-point>` - Destination the data disk mounted to in VM.

```hocon
 default-runtime-attributes {
   dataDisk: "<disk-type> <disk-size-in-GB> <mount-point>"
 }
```

###CallCaching
BCS supports CallCaching feature when the docker image is from Alibaba Cloud Container Registry.
The configuration file will look like the following:
```hocon
call-caching {
  enabled = true
  invalidate-bad-cache-results = true

}

docker {
  hash-lookup {
    enabled = true
    method = "remote"
    alibabacloudcr {
      num-threads = 5
      auth {
        access-id = xxxx
        access-key = yyyy
        security-token = zzzz
      }
    }
  }
}

backend {
  providers {
    BCS {
      config {
        # BCS and OSS related configurations mentioned above
        filesystems {
          oss {
            caching {
               duplication-strategy = "reference"
               invalidate-bad-cache-results = true
            }
            # ... to be filled in
          }
        }
        default-runtime-attributes {
          docker: "registry.cn-shanghai.aliyuncs.com/batchcompute/myubuntu:0.2"
          # ... to be filled in
        }
      }
    }
  }
}
```

- `docker.hash-lookup.method` - BCS only supports `remote` method for hash-lookup
- `filesystems.oss.caching.duplication-strategy` - BCS only supports `reference` for duplication strategy.