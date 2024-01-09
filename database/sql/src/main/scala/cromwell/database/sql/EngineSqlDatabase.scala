package cromwell.database.sql

trait EngineSqlDatabase
    extends SqlDatabase
    with JobKeyValueSqlDatabase
    with CallCachingSqlDatabase
    with JobStoreSqlDatabase
    with WorkflowStoreSqlDatabase
    with SubWorkflowStoreSqlDatabase
    with DockerHashStoreSqlDatabase
