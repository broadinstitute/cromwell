# AWS Batch backend (beta)

Check out the [getting started guide](../tutorials/AwsBatch101.md) for the bulk of our documentation.

AWS support is fairly new to Cromwell and this reference section will expand as features are added and documented.

### Disks

Cromwell performs automatic disk sizing on your behalf when running with the AWS backend, so attributes like 
```
disks: "local-disk"
```
or
```
disks: "/some/mnt"
```
are adequate to specify a disk that will suffice to complete your task.

To facilitate the running of workflows originally authored for Pipelines API on Google Cloud Platform, Cromwell's AWS backend can also interpret attributes like
```
disks: "local-disk 20 SSD"
```
and
```
disks: "/some/mnt 20 SSD"
```
The size information and HDD/SSD have no effect on this backend and Cromwell simply drops them.
