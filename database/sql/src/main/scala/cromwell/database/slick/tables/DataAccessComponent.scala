package cromwell.database.slick.tables

abstract class DataAccessComponent extends DriverComponent {
  def schema: driver.SchemaDescription
}
