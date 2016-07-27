package cromwell.database.sql

trait SqlDatabase extends AutoCloseable
  with OldeWorldeSqlDatabase
  with MetadataSqlDatabase
  with WorkflowStoreSqlDatabase
  with BackendKVStoreSqlDatabase
  with JobStoreSqlDatabase
  with BackendKVStoreSqlDatabase
