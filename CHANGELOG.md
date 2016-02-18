# Cromwell Change Log

## 0.18

* The deprecated parse, validate, inputs and highlight functionality from the command line tool has been removed in favor of wdltool (https://github.com/broadinstitute/wdltool) 
* Workflow options can now include a `defaultRuntimeOptions` section, so that the same runtime attribute is not needed in every single WDL task. E.g.:
```
{
  "defaultRuntimeOptions": {
    "docker": "ubuntu:latest"
  }
}
```
* Changed the JES runtime attributes `defaultDisks` and `defaultZones` to be simply `disks` and `zones` respectively.
* Liquibase scripts now run automatically. Non-persistent, in-memory databases are not affected. However Cromwell will
not start if it detects evidence of manually run liquibase migrations in a persistent database. Instead, before Cromwell
will start cleanly, the database should backed up, and then this SQL should be manually executed:
```sql
update DATABASECHANGELOG
  set MD5SUM = null,
    FILENAME = substr(FILENAME, instr(FILENAME, "src/main/migrations/") + length("src/main/migrations/"))
  where FILENAME like '%src/main/migrations/%'
```
* Added Preemptible VMs support for JES. This has impacts on the API Endpoint responses as a Call/Shard can now be attempted multiple times. Each attempt will have its own entry.
* Added custom thread pool to workaround Slick [deadlocks](https://github.com/slick/slick/issues/1274). The thread pool
size defaults to the Slick configuration value `db.numThreads`, but may be increased up to Slick's
`db.maxConnections`, via a new property `actionThreadPoolSize`.
* Added support for [size](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#float-sizefile-string) WDL standard library function.

* Allow for runtime attribute values to be interpreted as full expressions.  For example:
```
task example {
  String ubuntu_tag
  command { ... }
  runtime {
    docker: "ubuntu:" + ubuntu_tag
  }
}
```
* Add runtime attributes in Call metadata :
```
{
  "workflowName": "hello",
  "calls": {
    "hello.hello": [
      {
        ...,
        "runtimeAttributes": {
                  "preemptible": "0",
                  "failOnStderr": "false",
                  "disks": "local-disk 10 SSD",
                  "continueOnReturnCode": "0",
                  "docker": "ubuntu:latest",
                  "cpu": "1",
                  "zones": "us-central1-a",
                  "memory": "2GB"
                },
        ... 
      }
    ]
  }
}
```
* Added "preemptible" field in Call metadata. This only appears if the backend is JES.
```
{
  "workflowName": "hello",
  "calls": {
    "hello.hello": [
      {
        ...,
        "preemptible": "true"
        ... 
      }
    ]
  }
}
```
* Bug fix: Tasks that changed directory would fail on JES because their return code file was written to the new directory instead of an absolute path
* Bug fix: Using `write_*` functions in a Task's command (e.g. `./my_script --file=${write_file(my_array)}`) will now work with JES
* Changing format of the 'disks' runtime attribute slightly to allow for mounting disks at specific mountpoints
```
task disk_test {
  command { ... }
  runtime {
    disks: "local-disk 20 SSD, /mnt/mnt1 200 HDD"
  }
}
```
