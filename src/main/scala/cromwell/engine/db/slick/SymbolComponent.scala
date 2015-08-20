package cromwell.engine.db.slick

import java.sql.Clob

case class Symbol
(
  workflowExecutionId: Int,
  scope: String,
  name: String,
  index: Int, // https://bugs.mysql.com/bug.php?id=8173
  io: String,
  wdlType: String,
  wdlValue: Option[Clob],
  symbolId: Option[Int] = None
  )

trait SymbolComponent {
  this: DriverComponent with WorkflowExecutionComponent =>

  import driver.api._

  class Symbols(tag: Tag) extends Table[Symbol](tag, "SYMBOL") {
    def symbolId = column[Int]("SYMBOL_ID", O.PrimaryKey, O.AutoInc)

    def workflowExecutionId = column[Int]("WORKFLOW_EXECUTION_ID")

    def scope = column[String]("SCOPE")

    def name = column[String]("NAME")

    def index = column[Int]("INDEX")

    def io = column[String]("IO")

    def wdlType = column[String]("WDL_TYPE")

    def wdlValue = column[Option[Clob]]("WDL_VALUE")

    override def * = (workflowExecutionId, scope, name, index, io, wdlType, wdlValue, symbolId.?) <>
      (Symbol.tupled, Symbol.unapply)

    def workflowExecution = foreignKey(
      "FK_SYMBOL_WORKFLOW_EXECUTION_ID", workflowExecutionId, workflowExecutions)(_.workflowExecutionId)

    def uniqueKey = index("UK_SYM_WORKFLOW_EXECUTION_ID_SCOPE_NAME_ITERATION_IO",
      (workflowExecutionId, scope, name, index, io), unique = true)
  }

  protected val symbols = TableQuery[Symbols]

  val symbolsAutoInc = symbols returning symbols.
    map(_.symbolId) into ((a, id) => a.copy(symbolId = Some(id)))

  // Convenience function
  def symbolsByWorkflowExecutionUuidAndIoAndMaybeScope(workflowExecutionUuid: String,
                                                       io: String, scopeOption: Option[String]) = {
    scopeOption match {
      case Some(scope) => symbolsByWorkflowExecutionUuidAndIoAndScope(workflowExecutionUuid, io, scope)
      case None => symbolsByWorkflowExecutionUuidAndIo(workflowExecutionUuid, io)
    }
  }

  def symbolsByScopeAndName(workflowExecutionUuid: String, scope: String, name: String) = symbolsByUuidAndScopeAndName(workflowExecutionUuid, scope, name)

  val allSymbols = Compiled(
    (workflowExecutionUuid: Rep[String]) => for {
      symbol <- symbols
      workflowExecution <- symbol.workflowExecution
      if workflowExecution.workflowExecutionUuid === workflowExecutionUuid
    } yield symbol)

  val symbolsByUuidAndScopeAndName = Compiled(
    (workflowExecutionUuid: Rep[String], scope: Rep[String], name: Rep[String]) => for {
      symbol <- symbols
      if symbol.scope === scope && symbol.name === name
      workflowExecution <- symbol.workflowExecution
      if workflowExecution.workflowExecutionUuid === workflowExecutionUuid
    } yield symbol)

  val symbolsByWorkflowExecutionUuidAndIo = Compiled(
    (workflowExecutionUuid: Rep[String], io: Rep[String]) => for {
      symbol <- symbols
      if symbol.io === io
      workflowExecution <- symbol.workflowExecution
      if workflowExecution.workflowExecutionUuid === workflowExecutionUuid
    } yield symbol)

  val symbolsByWorkflowExecutionUuidAndIoAndScope = Compiled(
    (workflowExecutionUuid: Rep[String], io: Rep[String], scope: Rep[String]) => for {
      symbol <- symbols
      if symbol.io === io && symbol.scope === scope
      workflowExecution <- symbol.workflowExecution
      if workflowExecution.workflowExecutionUuid === workflowExecutionUuid
    } yield symbol)
}
