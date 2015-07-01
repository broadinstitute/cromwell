package cromwell.engine.db.slick

case class Symbol
(
  workflowExecutionId: Int,
  scope: String,
  name: String,
  iteration: Int, // https://bugs.mysql.com/bug.php?id=8173
  io: String,
  wdlType: String,
  wdlValue: Option[String],
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

    def iteration = column[Int]("ITERATION")

    def io = column[String]("IO")

    def wdlType = column[String]("WDL_TYPE")

    def wdlValue = column[Option[String]]("WDL_VALUE")

    override def * = (workflowExecutionId, scope, name, iteration, io, wdlType, wdlValue, symbolId.?) <>
      (Symbol.tupled, Symbol.unapply)

    def workflow = foreignKey(
      "FK_SYMBOL_WORKFLOW_EXECUTION_ID", workflowExecutionId, workflowExecutions)(_.workflowExecutionId)

    def uniqueKey = index("UK_SYM_WORKFLOW_EXECUTION_ID_SCOPE_NAME_ITERATION_IO",
      (workflowExecutionId, scope, name, iteration, io), unique = true)
  }

  val symbols = TableQuery[Symbols]

  val symbolsAutoInc = symbols returning symbols.
    map(_.symbolId) into ((a, id) => a.copy(symbolId = Some(id)))

  val symbolByID = Compiled(
    (id: Rep[Int]) => for {
      symbol <- symbols
      if symbol.symbolId === id
    } yield symbol)
}
