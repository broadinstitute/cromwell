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

