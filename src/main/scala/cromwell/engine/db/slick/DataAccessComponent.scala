package cromwell.engine.db.slick

import slick.driver.JdbcProfile

import scala.reflect.runtime._

class DataAccessComponent(val driver: JdbcProfile)
  extends DriverComponent
  with ExecutionComponent
  with JesJobComponent
  with LocalJobComponent
  with SymbolComponent
  with WorkflowExecutionComponent {

  import driver.api._

  def this(driverName: String) {
    this(DataAccessComponent.getObject[JdbcProfile](driverName))
  }

  def insertWorkflowExecution(workflowExecution: WorkflowExecution): DBIO[WorkflowExecution] = {
    workflowExecutionsAutoInc += workflowExecution
  }

  def insertSymbols(symbols: Seq[Symbol]): DBIO[Seq[Symbol]] = {
    symbolsAutoInc ++= symbols
  }
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
