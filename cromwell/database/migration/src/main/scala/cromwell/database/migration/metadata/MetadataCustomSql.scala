package cromwell.database.migration.metadata

import java.time.OffsetDateTime

import liquibase.change.custom.CustomSqlChange
import liquibase.database.Database
import liquibase.exception.ValidationErrors
import liquibase.resource.ResourceAccessor
import liquibase.statement.SqlStatement
import liquibase.statement.core.RawSqlStatement

object MetadataCustomSql {
  val Offset = OffsetDateTime.now().getOffset.toString
}

abstract class MetadataCustomSql extends CustomSqlChange {

  def queries: Array[String]

  override def generateStatements(database: Database): Array[SqlStatement] = {
    queries map { query =>  new RawSqlStatement(query) }
  }

  override def setUp(): Unit = ()

  override def validate(database: Database): ValidationErrors = {
    new ValidationErrors()
  }

  override def setFileOpener(resourceAccessor: ResourceAccessor): Unit = ()
}
