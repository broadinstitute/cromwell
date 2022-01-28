# HTTP Inputs

## Overview

For shared filesystem and Google Pipelines API (PAPI) version 2 backends Cromwell can support workflow inputs specified by `http` and `https` URLs.
Please note this is not true "filesystem" support for HTTP URLs;
if inputs to a workflow are specified by HTTP URLs the outputs of steps will nevertheless appear at local or GCS paths and not HTTP
URLs.

### Configuration

Cromwell's default configuration defines an instance of the HTTP filesystem named `http`. There is no additional configuration
required for the HTTP filesystem itself so adding HTTP filesystem support to a backend is a simple as
adding a reference to this filesystem within the backend's `filesystems` stanza. e.g. Cromwell's default `Local` shared filesystem
backend is configured like this (a PAPI version 2 backend would be configured in a similar way):

```
backend {
  default = "Local"
  providers {
    Local {
      ...
      config {
        filesystems {
          local {
            ...
          }
          http { }
        }
      }
      ...
    }
    ...
  }
}
```

If there is a need to turn off this `http` filesystem in the default `Local` backend the following Java property
allows for this: `-Dbackend.providers.Local.config.filesystems.http.enabled=false`.

### Caveats

Using HTTP inputs in Cromwell can produce some unexpected behavior:
- Files specified by HTTP URIs will be renamed locally, so programs that rely on file extensions or other filenaming conventions may not function properly.
- Files located in the same remote HTTP-defined directory will not be colocated locally. This can cause problems if a program is expecting an index file (e.g. `.fai`) to appear in the same directory as the associated data file (e.g. `.fa`) without specifying the index location.
