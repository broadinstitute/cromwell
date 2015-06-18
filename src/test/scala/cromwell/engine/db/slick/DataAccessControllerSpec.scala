package cromwell.engine.db.slick

import java.sql.Timestamp
import java.util.UUID

import cromwell.engine.store.ExecutionStore.ExecutionStatus
import org.scalatest.concurrent.ScalaFutures
import org.scalatest.{FlatSpec, Matchers}

import scala.concurrent.Future

class DataAccessControllerSpec extends FlatSpec with Matchers with ScalaFutures {
  SlickTestDatabase.checkInitialized()
  import DataAccessController.dataAccess.driver.api._

  "DataAccessController" should "access the database" in {
    DataAccessController.database.close()
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
