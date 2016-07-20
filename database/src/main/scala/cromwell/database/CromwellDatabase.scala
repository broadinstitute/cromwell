package cromwell.database

import cromwell.database.slick.SlickDatabase

object CromwellDatabase {
  val databaseInterface: SqlDatabase = new SlickDatabase()
}

trait Database {
  def databaseInterface: SqlDatabase
}

trait CromwellDatabase extends Database {
  override def databaseInterface: SqlDatabase = CromwellDatabase.databaseInterface
}
