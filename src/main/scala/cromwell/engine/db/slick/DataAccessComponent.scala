package cromwell.engine.db.slick

import slick.driver.JdbcProfile

import scala.reflect.runtime._

class DataAccessComponent(val driver: JdbcProfile)
  extends DriverComponent
  with WorkflowExecutionComponent
  with WorkflowExecutionAuxComponent
  with SymbolComponent
  with ExecutionComponent
  with JesJobComponent
  with LocalJobComponent {

  import driver.api._

  def this(driverName: String) {
    this(DataAccessComponent.getObject[JdbcProfile](driverName))
  }

  lazy val schema =
    workflowExecutions.schema ++
      workflowExecutionAuxes.schema ++
      symbols.schema ++
      executions.schema ++
      localJobs.schema ++
      jesJobs.schema
}

object DataAccessComponent {
  // TODO: move to DSDE common util?
  private def getObject[T](objectName: String): T = {
    // via
    //   http://stackoverflow.com/questions/23466782/scala-object-get-reference-from-string-in-scala-2-10
    //   https://github.com/anvie/slick-test/blob/045f4db610d3b91bf928a53f2bc7b6ae17c35985/slick-util/src/main/scala/scala/slick/codegen/ModelGenerator.scala
    val staticModule = currentMirror.staticModule(objectName)
    val reflectModule = currentMirror.reflectModule(staticModule)
    val instance = reflectModule.instance
    instance.asInstanceOf[T]
  }
}
