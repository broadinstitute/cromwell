package cromwell.engine.db.slick

import java.sql.Timestamp
import java.util.{Date, UUID}

import cromwell.binding.command.Command
import cromwell.binding.types.WdlStringType
import cromwell.binding.values.WdlString
import cromwell.binding.{Call, Task}
import cromwell.engine.db.LocalCallInfo
import cromwell.engine.store.ExecutionStore.ExecutionStatus
import cromwell.engine.store.SymbolStore.SymbolStoreKey
import cromwell.engine.store.SymbolStoreEntry
import cromwell.engine.{WorkflowRunning, WorkflowSubmitted}
import cromwell.util.SampleWdl.HelloWorld
import org.scalatest.concurrent.ScalaFutures
import org.scalatest.{FlatSpec, Matchers}

import scala.concurrent.Future

class DataAccessControllerSpec extends FlatSpec with Matchers with ScalaFutures {
  SlickTestDatabase.checkInitialized()

  import DataAccessController.dataAccess.driver.api._

  "DataAccessController" should "access the database" in {
    DataAccessController.database.close()
  }

  val workflowIdForCreate = UUID.randomUUID()

  it should "create a workflow for just reading" in {
    val key = new SymbolStoreKey("myScope", "myName", None, true)
    val entries = Seq(new SymbolStoreEntry(key, WdlStringType, Option(new WdlString("testStringValue"))))
    DataAccessController.createWorkflow(workflowIdForCreate, "createUri", HelloWorld.WdlSource, HelloWorld.JsonInputs, entries)
  }

  it should "retrieve the workflow for just reading" in {
    val results = DataAccessController.query(workflowId = Some(Seq(workflowIdForCreate)))
    results should have size 1

    val workflowResult = results.head
    workflowResult.state should be(WorkflowSubmitted)
    workflowResult.wdlUri should be("createUri")
    workflowResult.startTime.getTime should be(new Date().getTime +- (60 * 1000L))
    workflowResult.endTime should be(empty)
    workflowResult.wdlSource shouldEqual HelloWorld.WdlSource
    workflowResult.jsonInputs shouldEqual HelloWorld.JsonInputs

    val resultCalls = workflowResult.calls
    resultCalls shouldBe empty

    val resultSymbols = workflowResult.symbols
    resultSymbols.size should be(1)
    val resultSymbol = resultSymbols.head
    val resultSymbolStoreKey = resultSymbol.key
    resultSymbolStoreKey.scope should be("myScope")
    resultSymbolStoreKey.name should be("myName")
    resultSymbolStoreKey.iteration should be(None)
    resultSymbolStoreKey.input should be(right = true) // IntelliJ highlighting
    resultSymbol.wdlType should be(WdlStringType)
    resultSymbol.wdlValue shouldNot be(empty)
    resultSymbol.wdlValue.get should be(new WdlString("testStringValue"))
  }

  it should "update a workflow state" in {
    val workflowId = UUID.randomUUID()
    val key = new SymbolStoreKey("myScope", "myName", None, true)
    val entries = Seq(new SymbolStoreEntry(key, WdlStringType, Option(new WdlString("testStringValue"))))
    DataAccessController.createWorkflow(workflowId, "createUri", HelloWorld.WdlSource, HelloWorld.JsonInputs, entries)
    DataAccessController.updateWorkflow(workflowId, WorkflowRunning)
    val results = DataAccessController.query(workflowId = Some(Seq(workflowId)))
    results.size should be(1)
    results.head.state should be(WorkflowRunning)
  }

  it should "add and update a workflow call" in {
    val workflowId = UUID.randomUUID()
    val key = new SymbolStoreKey("myScope", "myName", None, true)
    val entries = Seq(new SymbolStoreEntry(key, WdlStringType, Option(new WdlString("testStringValue"))))
    val task = new Task("taskName", new Command(Seq.empty), Seq.empty, Map.empty)
    val call = new Call(None, "fully.qualified.name", task, Map.empty, null)

    DataAccessController.createWorkflow(workflowId, "createUri", HelloWorld.WdlSource, HelloWorld.JsonInputs, entries)
    DataAccessController.updateWorkflow(workflowId, WorkflowRunning)

    val callInfoInsert = new LocalCallInfo(call.fullyQualifiedName, ExecutionStatus.NotStarted, -1, "test command", -1)
    DataAccessController.updateCall(workflowId, call, None, Option(callInfoInsert), None)
    val insertResults = DataAccessController.query(workflowId = Some(Seq(workflowId)))
    insertResults.size should be(1)
    insertResults.head.state should be(WorkflowRunning)
    val insertResultCalls = insertResults.head.calls
    insertResultCalls.size should be(1)
    val insertResultCall = insertResultCalls.toSeq.head
    insertResultCall should be (a [LocalCallInfo])
    val insertResultLocalCall = insertResultCall.asInstanceOf[LocalCallInfo]
    insertResultLocalCall.callFqn should be("fully.qualified.name")
    insertResultLocalCall.status should be(ExecutionStatus.NotStarted)
    insertResultLocalCall.processId should be(-1)
    insertResultLocalCall.command should be("test command")
    insertResultLocalCall.resultCode should be(-1)

    val callInfoUpdate = new LocalCallInfo(call.fullyQualifiedName, ExecutionStatus.Running, 123, "test updated", 1234)
    DataAccessController.updateCall(workflowId, call, None, Option(callInfoUpdate), None)
    val updateResults = DataAccessController.query(workflowId = Some(Seq(workflowId)))
    updateResults.size should be(1)
    updateResults.head.state should be(WorkflowRunning)
    val updateResultCalls = updateResults.head.calls
    updateResultCalls.size should be(1)
    val updateResultCall = updateResultCalls.toSeq.head
    updateResultCall should be (a [LocalCallInfo])
    val updateResultLocalCall = updateResultCall.asInstanceOf[LocalCallInfo]
    updateResultLocalCall.callFqn should be("fully.qualified.name")
    updateResultLocalCall.status should be(ExecutionStatus.Running)
    updateResultLocalCall.processId should be(123)
    updateResultLocalCall.command should be("test updated")
    updateResultLocalCall.resultCode should be(1234)
  }

  var workflowExecutionId: Option[Int] = None
  var symbolId: Option[Int] = None
  val uuid: String = UUID.randomUUID.toString

  it should "create a workflow execution" in {
    val insert = DataAccessController.dataAccess.workflowExecutionsAutoInc +=
      new WorkflowExecution(uuid, "wdl_uri", ExecutionStatus.Starting.toString, new Timestamp(1000))

    val workflowExecutionFuture: Future[WorkflowExecution] = DataAccessController.database.run(insert.transactionally)

    whenReady(workflowExecutionFuture) { workflowExecution =>
      workflowExecution.workflowExecutionId shouldNot be(empty)
      workflowExecution.workflowExecutionUuid should be(uuid)
      workflowExecution.wdlUri should be("wdl_uri")
      workflowExecution.status should be(ExecutionStatus.Starting.toString)
      workflowExecution.startDt should be(new Timestamp(1000))
      workflowExecution.endDt should be(empty)
      workflowExecutionId = workflowExecution.workflowExecutionId
    }
  }

  it should "read a workflow execution" in {
    assume(workflowExecutionId.isDefined)

    val select = DataAccessController.dataAccess.workflowExecutions.filter(
      _.workflowExecutionId === workflowExecutionId).result.head
    val workflowExecutionFuture: Future[WorkflowExecution] = DataAccessController.database.run(select.transactionally)

    whenReady(workflowExecutionFuture) { workflowExecution =>
      workflowExecution.workflowExecutionId should be(workflowExecutionId)
      workflowExecution.workflowExecutionUuid should be(uuid)
      workflowExecution.wdlUri should be("wdl_uri")
      workflowExecution.status should be(ExecutionStatus.Starting.toString)
      workflowExecution.startDt should be(new Timestamp(1000))
      workflowExecution.endDt should be(empty)
    }
  }

  it should "create a symbol" in {
    assume(workflowExecutionId.isDefined)

    val insert = DataAccessController.dataAccess.symbolsAutoInc +=
      new Symbol(workflowExecutionId.get, "symbol_scope", "symbol_name", None, "OUTPUT", "String", "hello world")

    val symbolFuture: Future[Symbol] = DataAccessController.database.run(insert.transactionally)

    whenReady(symbolFuture) { symbol =>
      symbol.symbolId shouldNot be(empty)
      symbol.scope should be("symbol_scope")
      symbol.name should be("symbol_name")
      symbol.iteration should be(empty)
      symbol.io should be("OUTPUT")
      symbol.wdlType should be("String")
      symbol.wdlValue should be("hello world")
      symbolId = symbol.symbolId
    }
  }

  it should "read a symbol" in {
    assume(symbolId.isDefined)
    val select = DataAccessController.dataAccess.symbols.filter(
      _.symbolId === symbolId).result.head
    val symbolFuture: Future[Symbol] = DataAccessController.database.run(select.transactionally)

    whenReady(symbolFuture) { symbol =>
      symbol.symbolId should be(symbolId)
      symbol.scope should be("symbol_scope")
      symbol.name should be("symbol_name")
      symbol.iteration should be(empty)
      symbol.io should be("OUTPUT")
      symbol.wdlType should be("String")
      symbol.wdlValue should be("hello world")
    }
  }
}
