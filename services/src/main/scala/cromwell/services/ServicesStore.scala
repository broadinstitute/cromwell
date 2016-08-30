package cromwell.database

import cromwell.database.slick.SlickDatabase
import cromwell.database.sql.SqlDatabase


object CromwellDatabase {
  val databaseInterface: SqlDatabase = new SlickDatabase()
}

trait Database {
  def databaseInterface: SqlDatabase
}

trait CromwellDatabase extends Database {
  override def databaseInterface: SqlDatabase = CromwellDatabase.databaseInterface
}
