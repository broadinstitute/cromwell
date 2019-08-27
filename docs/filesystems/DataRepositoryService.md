# Data Repository Service (DRS)

The Cromwell configuration for DRS is as follows:

**Filesystem Configuration**

```hocon
drs {
    # A reference to a potentially different auth required to contact Martha service.
    auth = "application-default"
}
```

The `auth` field refers to the authentication schema that should be used to authenticate requests to Martha service.

The `drs` section needs to be added to
- `engine.filesystems` block
- [PapiV2](http://cromwell.readthedocs.io/en/develop/backends/Google) backend's `filesystems` block


**Localization Configuration**

DRS localization must be configured with the docker image to use.

```hocon
drs {
    localization {
        docker-image = "broadinstitute/drs-localizer:latest"
    }
}
```


**Example**

A sample configuration for DRS filesystem might look like:

```hocon
engine {
  filesystems {
    # ... other filesystems here, probably gcs, and then ...
    drs {
      auth = "application-default"
    }
  }
}

backend {
    # ... other global backend config here, probably just setting the default ...
    providers {
        # ... other providers here ...
        Papi {
            actor-factory = "cromwell.backend.google.pipelines.v2alpha1.PipelinesApiLifecycleActorFactory"
            config {
                # ... other config here ...
                filesystems {
                    # ... other filesystems here, probably gcs, and then ...
                    drs {
                        auth = "application-default"
                    }
                }
            }
        }
    }
}

drs {
    localization {
        docker-image = "broadinstitute/drs-localizer:latest"
    }
}
```
