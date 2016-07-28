package cromwell.database.migration.metadata.table.symbol

import java.sql.ResultSet

class ResultSetIterator(rs: ResultSet) extends Iterator[ResultSet] {
  def hasNext: Boolean = rs.next()
  def next(): ResultSet = rs
}
