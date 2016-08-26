package cromwell.database.sql.tables

trait DatabaseSimpleton {
  def simpletonKey: String
  def simpletonValue: String
  def wdlType: String
}
