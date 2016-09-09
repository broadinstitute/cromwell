package cromwell.database.migration.metadata.table.symbol

import java.sql.{PreparedStatement, ResultSet}

import liquibase.database.jvm.JdbcConnection

class QueryPaginator(statement: PreparedStatement,
                     batchSize: Int,
                     count: Int) extends Iterator[ResultSet] {
  var cursor = 0

  def next(): ResultSet =  {
    statement.setInt(1, cursor)
    statement.setInt(2, cursor + batchSize)

    cursor += batchSize
    statement.executeQuery()
  }

  def hasNext(): Boolean = cursor <= count
}
