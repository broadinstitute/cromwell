# Cromwell Change Log

## 0.18

* Workflow options can now include a `defaultRuntimeOptions` section, so that the same runtime attribute is not needed in every single WDL task. E.g.:
```
{
  "defaultRuntimeOptions": {
    "docker": "ubuntu:latest"
  }
}
```
* Changed the JES runtime attributes `defaultDisks` and `defaultZones` to be simply `disks` and `zones` respectively.

* Added Preemptible VMs support for JES. This has impacts on the API Endpoint responses as a Call/Shard can now be attempted multiple times. Each attempt will have its own entry.

