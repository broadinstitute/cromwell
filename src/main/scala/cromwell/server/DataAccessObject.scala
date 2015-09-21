package cromwell.server

import cromwell.engine.db.DataAccess

trait DataAccessObject {
  lazy val dataAccess: DataAccess = DataAccess()
}
