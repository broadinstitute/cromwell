package cromwell.database.slick

import java.sql.Clob

import cromwell.database.obj.Symbol

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

    def reportableResult = column[Boolean]("REPORTABLE_RESULT")

    def wdlType = column[String]("WDL_TYPE")

    def wdlValue = column[Option[Clob]]("WDL_VALUE")

    def symbolHash = column[Option[String]]("HASH")

    override def * = (workflowExecutionId, scope, name, index, io, reportableResult, wdlType, wdlValue, symbolHash, symbolId.?) <>
      (Symbol.tupled, Symbol.unapply)

    def workflowExecution = foreignKey(
      "FK_SYMBOL_WORKFLOW_EXECUTION_ID", workflowExecutionId, workflowExecutions)(_.workflowExecutionId)

    def uniqueKey = index("UK_SYM_WORKFLOW_EXECUTION_ID_SCOPE_NAME_ITERATION_IO",
      (workflowExecutionId, scope, name, index, io), unique = true)
  }

  protected val symbols = TableQuery[Symbols]

  val symbolIdsAutoInc = symbols returning symbols.map(_.symbolId)

  val symbolsByWorkflowExecutionUuid = Compiled(
    (workflowExecutionUuid: Rep[String]) => for {
      symbol <- symbols
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

  val symbolsByWorkflowExecutionUuidAndIoAndScopeAndIndex = Compiled(
    (workflowExecutionUuid: Rep[String], io: Rep[String], scope: Rep[String], index: Rep[Int]) => for {
      symbol <- symbols
      if symbol.io === io
      if symbol.scope === scope
      if symbol.index === index
      workflowExecution <- symbol.workflowExecution
      if workflowExecution.workflowExecutionUuid === workflowExecutionUuid
    } yield symbol)

  val symbolWdlTypeAndWdlValueByWorkflowAndScopeAndIndexAndName = Compiled(
    (workflowExecutionId: Rep[Int], scope: Rep[String], index: Rep[Int], name: Rep[String]) => for {
      symbol <- symbols
      if symbol.workflowExecutionId === workflowExecutionId
      if symbol.scope === scope
      if symbol.index === index
      if symbol.name === name
    } yield (symbol.wdlType, symbol.wdlValue))

  val symbolsForWorkflowOutput = Compiled(
    (workflowExecutionUuid: Rep[String]) => for {
      symbol <- symbols
      if symbol.reportableResult === true
      workflowExecution <- symbol.workflowExecution
      if workflowExecution.workflowExecutionUuid === workflowExecutionUuid
    } yield symbol)
}
