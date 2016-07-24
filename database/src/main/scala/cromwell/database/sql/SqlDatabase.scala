package cromwell.database.sql

trait SqlDatabase extends AutoCloseable
  with OldeWorldeSqlDatabase
  with MetadataSqlDatabase
  with WorkflowStoreSqlDatabase
