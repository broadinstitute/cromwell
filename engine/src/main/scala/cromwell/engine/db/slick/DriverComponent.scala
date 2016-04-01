package cromwell.engine.db.slick

import slick.driver.JdbcProfile

trait DriverComponent {
  val driver: JdbcProfile
}
