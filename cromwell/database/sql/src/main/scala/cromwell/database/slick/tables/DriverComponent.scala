package cromwell.database.slick.tables

import slick.jdbc.JdbcProfile

trait DriverComponent {
  val driver: JdbcProfile
}
