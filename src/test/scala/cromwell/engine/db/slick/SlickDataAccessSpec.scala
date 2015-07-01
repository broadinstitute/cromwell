package cromwell.engine.db.slick

import java.sql.SQLException
import java.util.UUID

import cromwell.binding._
import cromwell.binding.command.Command
import cromwell.binding.types.WdlStringType
import cromwell.binding.values.WdlString
import cromwell.engine._
import cromwell.engine.backend.local.LocalBackend
import cromwell.engine.db.DataAccess.WorkflowInfo
import cromwell.engine.db.LocalCallBackendInfo
import org.scalatest.concurrent.ScalaFutures
import org.scalatest.time.{Millis, Seconds, Span}
import org.scalatest.{FlatSpec, Matchers}

import scala.concurrent.ExecutionContext

class SlickDataAccessSpec extends FlatSpec with Matchers with ScalaFutures {
  implicit val ec = ExecutionContext.global

  implicit val defaultPatience = PatienceConfig(timeout = Span(5, Seconds), interval = Span(100, Millis))

  // Tests against main database used for command line
  "SlickDataAccess (main.hsqldb)" should behave like databaseWithConfig("main.hsqldb", testRequired = true)

  // Tests using liquibase, but in memory
  "SlickDataAccess (test.hsqldb)" should behave like databaseWithConfig("test.hsqldb", testRequired = true)

  // If able to connect, then also run the tests on mysql, but it's not required
  "SlickDataAccess (test.mysql)" should behave like databaseWithConfig("test.mysql", testRequired = false)

  def databaseWithConfig(path: => String, testRequired: => Boolean): Unit = {

    lazy val testDatabase = new TestSlickDatabase(path)
    lazy val dataAccess = testDatabase.slickDataAccess
    lazy val canConnect = {
      testRequired || (DatabaseConfig.rootDatabaseConfig.hasPath(path) && testDatabase.isValidConnection.futureValue)
    }

    it should "setup via liquibase if necessary" in {
      assume(canConnect || testRequired)
      if (testDatabase.useLiquibase)
        testDatabase.setupLiquibase()
    }

    it should "create and retrieve the workflow for just reading" in {
      assume(canConnect || testRequired)
      val workflowId = UUID.randomUUID()
      val workflowInfo = new WorkflowInfo(workflowId, "source", "{}")

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, Seq.empty, Seq.empty, new LocalBackend)
        _ <- dataAccess.getWorkflowsByState(Seq(WorkflowSubmitted)) map { results =>
          results shouldNot be(empty)

          val workflowResultOption = results.find(_.workflowId == workflowId)
          workflowResultOption shouldNot be(empty)
          val workflowResult = workflowResultOption.get
          workflowResult.workflowId should be(workflowId)
          workflowResult.wdlSource should be("source")
          workflowResult.wdlJson should be("{}")
        }
      } yield ()).futureValue
    }

    it should "fail when saving a workflow twice" in {
      assume(canConnect || testRequired)
      val workflowId = UUID.randomUUID()
      val workflowInfo = new WorkflowInfo(workflowId, "source", "{}")

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, Seq.empty, Seq.empty, new LocalBackend)
        _ <- dataAccess.createWorkflow(workflowInfo, Seq.empty, Seq.empty, new LocalBackend).failed map { ex =>
          ex should be(a[SQLException])
        }
      } yield ()).futureValue
    }

    it should "update and get a workflow state" in {
      assume(canConnect || testRequired)
      val workflowId = UUID.randomUUID()
      val workflowInfo = new WorkflowInfo(workflowId, "source", "{}")

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, Seq.empty, Seq.empty, new LocalBackend)
        _ <- dataAccess.updateWorkflowState(workflowId, WorkflowRunning)
        _ <- dataAccess.getWorkflowsByState(Seq(WorkflowRunning)) map { results =>
          results shouldNot be(empty)

          val workflowResultOption = results.find(_.workflowId == workflowId)
          workflowResultOption shouldNot be(empty)
          val workflowResult = workflowResultOption.get
          workflowResult.workflowId should be(workflowId)
          workflowResult.wdlSource should be("source")
          workflowResult.wdlJson should be("{}")
        }
      } yield ()).futureValue
    }

    it should "get a symbol input" in {
      assume(canConnect || testRequired)
      val callFqn = "call.fully.qualified.scope"
      val symbolFqn = "symbol.fully.qualified.scope"
      val workflowId = UUID.randomUUID()
      val workflowInfo = new WorkflowInfo(workflowId, "source", "{}")
      val key = new SymbolStoreKey(callFqn, symbolFqn, None, input = true)
      val entry = new SymbolStoreEntry(key, WdlStringType, Option(new WdlString("testStringValue")))
      val task = new Task("taskName", new Command(Seq.empty), Seq.empty, Map.empty)
      val call = new Call(None, callFqn, task, Map.empty, null)

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, Seq(entry), Seq.empty, new LocalBackend)
        _ <- dataAccess.updateWorkflowState(workflowId, WorkflowRunning)
        _ <- dataAccess.getInputs(workflowId, call) map { resultSymbols =>
          resultSymbols.size should be(1)
          val resultSymbol = resultSymbols.head
          val resultSymbolStoreKey = resultSymbol.key
          resultSymbolStoreKey.scope should be("call.fully.qualified.scope")
          resultSymbolStoreKey.name should be("symbol.fully.qualified.scope")
          resultSymbolStoreKey.iteration should be(None)
          resultSymbolStoreKey.input should be(right = true) // Inteillj highlighting
          resultSymbol.wdlType should be(WdlStringType)
          resultSymbol.wdlValue shouldNot be(empty)
          resultSymbol.wdlValue.get should be(new WdlString("testStringValue"))
        }
      } yield ()).futureValue
    }

    it should "get workflow state" in {
      assume(canConnect || testRequired)
      val workflowId = UUID.randomUUID()
      val workflowInfo = new WorkflowInfo(workflowId, "source", "{}")
      val key = new SymbolStoreKey("myScope", "myName", None, input = true)
      val entry = new SymbolStoreEntry(key, WdlStringType, Option(new WdlString("testStringValue")))

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, Seq(entry), Seq.empty, new LocalBackend)
        _ <- dataAccess.updateWorkflowState(workflowId, WorkflowSucceeded)
        _ <- dataAccess.getWorkflowState(workflowId) map { result =>
          result shouldNot be(empty)
          result.get should be(WorkflowSucceeded)
        }
      } yield ()).futureValue
    }

    it should "get initial execution status" in {
      assume(canConnect || testRequired)
      val workflowId = UUID.randomUUID()
      val workflowInfo = new WorkflowInfo(workflowId, "source", "{}")
      val task = new Task("taskName", new Command(Seq.empty), Seq.empty, Map.empty)
      val call = new Call(None, "fully.qualified.name", task, Map.empty, null)

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, Seq.empty, Seq(call), new LocalBackend)
        _ <- dataAccess.updateWorkflowState(workflowId, WorkflowRunning)
        _ <- dataAccess.getExecutionStatuses(workflowId) map { result =>
          result.size should be(1)
          val (fqn, status) = result.head
          fqn should be("fully.qualified.name")
          status should be(ExecutionStatus.NotStarted)
        }
      } yield ()).futureValue
    }

    it should "get updated execution status" in {
      assume(canConnect || testRequired)
      val workflowId = UUID.randomUUID()
      val workflowInfo = new WorkflowInfo(workflowId, "source", "{}")
      val task = new Task("taskName", new Command(Seq.empty), Seq.empty, Map.empty)
      val call = new Call(None, "fully.qualified.name", task, Map.empty, null)

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, Seq.empty, Seq(call), new LocalBackend)
        _ <- dataAccess.updateWorkflowState(workflowId, WorkflowRunning)
        _ <- dataAccess.setStatus(workflowId, Seq(call), ExecutionStatus.Running)
        _ <- dataAccess.getExecutionStatuses(workflowId) map { result =>
          result.size should be(1)
          val (fqn, status) = result.head
          fqn should be("fully.qualified.name")
          status should be(ExecutionStatus.Running)
        }
      } yield ()).futureValue
    }

    it should "fail to get an non-existent execution status" in {
      assume(canConnect || testRequired)
      val workflowIdNotSaved = UUID.randomUUID()

      (for {
        _ <- dataAccess.getExecutionStatuses(workflowIdNotSaved).failed map { ex =>
          ex should be(a[NoSuchElementException])
        }
      } yield ()).futureValue
    }

    it should "set and get an output" in {
      assume(canConnect || testRequired)
      val callFqn = "call.fully.qualified.scope"
      val symbolLqn = "symbol"
      val workflowId = UUID.randomUUID()
      val workflowInfo = new WorkflowInfo(workflowId, "source", "{}")
      val task = new Task("taskName", new Command(Seq.empty), Seq.empty, Map.empty)
      val call = new Call(None, callFqn, task, Map.empty, null)

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, Seq.empty, Seq.empty, new LocalBackend)
        _ <- dataAccess.updateWorkflowState(workflowId, WorkflowRunning)
        _ <- dataAccess.setOutputs(workflowId, call, Map(symbolLqn -> new WdlString("testStringValue")))
        _ <- dataAccess.getOutputs(workflowId) map { results =>
          results.size should be(1)
          val resultSymbol = results.head
          val resultSymbolStoreKey = resultSymbol.key
          resultSymbolStoreKey.scope should be("call.fully.qualified.scope")
          resultSymbolStoreKey.name should be("call.fully.qualified.scope.symbol")
          resultSymbolStoreKey.iteration should be(None)
          resultSymbolStoreKey.input should be(right = false) // Inteillj highlighting
          resultSymbol.wdlType should be(WdlStringType)
          resultSymbol.wdlValue shouldNot be(empty)
          resultSymbol.wdlValue.get should be(new WdlString("testStringValue"))
        }
      } yield ()).futureValue
    }

    it should "set and get an output by call" in {
      assume(canConnect || testRequired)
      val callFqn = "call.fully.qualified.scope"
      val symbolLqn = "symbol"
      val workflowId = UUID.randomUUID()
      val workflowInfo = new WorkflowInfo(workflowId, "source", "{}")
      val task = new Task("taskName", new Command(Seq.empty), Seq.empty, Map.empty)
      val call = new Call(None, callFqn, task, Map.empty, null)

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, Seq.empty, Seq.empty, new LocalBackend)
        _ <- dataAccess.updateWorkflowState(workflowId, WorkflowRunning)
        _ <- dataAccess.setOutputs(workflowId, call, Map(symbolLqn -> new WdlString("testStringValue")))
        _ <- dataAccess.getOutputs(workflowId, call) map { results =>
          results.size should be(1)
          val resultSymbol = results.head
          val resultSymbolStoreKey = resultSymbol.key
          resultSymbolStoreKey.scope should be("call.fully.qualified.scope")
          resultSymbolStoreKey.name should be("call.fully.qualified.scope.symbol")
          resultSymbolStoreKey.iteration should be(None)
          resultSymbolStoreKey.input should be(right = false) // Inteillj highlighting
          resultSymbol.wdlType should be(WdlStringType)
          resultSymbol.wdlValue shouldNot be(empty)
          resultSymbol.wdlValue.get should be(new WdlString("testStringValue"))
        }
      } yield ()).futureValue
    }

    it should "set and get the same symbol with IO as input then output" in {
      assume(canConnect || testRequired)
      val callFqn = "call.fully.qualified.scope"
      val symbolLqn = "symbol"
      val symbolFqn = callFqn + "." + symbolLqn
      val workflowId = UUID.randomUUID()
      val workflowInfo = new WorkflowInfo(workflowId, "source", "{}")
      val key = new SymbolStoreKey(callFqn, symbolFqn, None, input = true)
      val entry = new SymbolStoreEntry(key, WdlStringType, Option(new WdlString("testStringValue")))
      val task = new Task("taskName", new Command(Seq.empty), Seq.empty, Map.empty)
      val call = new Call(None, callFqn, task, Map.empty, null)

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, Seq(entry), Seq.empty, new LocalBackend)
        _ <- dataAccess.updateWorkflowState(workflowId, WorkflowRunning)
        _ <- dataAccess.setOutputs(workflowId, call, Map(symbolLqn -> new WdlString("testStringValue")))
        _ <- dataAccess.getOutputs(workflowId, call) map { results =>
          results.size should be(1)
          val resultSymbol = results.head
          val resultSymbolStoreKey = resultSymbol.key
          resultSymbolStoreKey.scope should be("call.fully.qualified.scope")
          resultSymbolStoreKey.name should be("call.fully.qualified.scope.symbol")
          resultSymbolStoreKey.iteration should be(None)
          resultSymbolStoreKey.input should be(right = false) // Inteillj highlighting
          resultSymbol.wdlType should be(WdlStringType)
          resultSymbol.wdlValue shouldNot be(empty)
          resultSymbol.wdlValue.get should be(new WdlString("testStringValue"))
        }
      } yield ()).futureValue
    }

    it should "fail when setting an existing symbol output" in {
      assume(canConnect || testRequired)
      val callFqn = "call.fully.qualified.scope"
      val symbolLqn = "symbol"
      val symbolFqn = callFqn + "." + symbolLqn
      val workflowId = UUID.randomUUID()
      val workflowInfo = new WorkflowInfo(workflowId, "source", "{}")
      val key = new SymbolStoreKey(callFqn, symbolFqn, None, input = false)
      val entry = new SymbolStoreEntry(key, WdlStringType, Option(new WdlString("testStringValue")))
      val task = new Task("taskName", new Command(Seq.empty), Seq.empty, Map.empty)
      val call = new Call(None, callFqn, task, Map.empty, null)

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, Seq(entry), Seq.empty, new LocalBackend)
        _ <- dataAccess.updateWorkflowState(workflowId, WorkflowRunning)
        _ <- dataAccess.setOutputs(workflowId, call, Map(symbolLqn -> new WdlString("testStringValue"))).failed map {
          ex =>
            ex should be(a[SQLException])
        }
      } yield ()).futureValue
    }

    it should "set and get a backend info" in {
      assume(canConnect || testRequired)
      val workflowId = UUID.randomUUID()
      val workflowInfo = new WorkflowInfo(workflowId, "source", "{}")
      val task = new Task("taskName", new Command(Seq.empty), Seq.empty, Map.empty)
      val call = new Call(None, "fully.qualified.name", task, Map.empty, null)
      val backendInfo = new LocalCallBackendInfo(ExecutionStatus.Running, Option(123), Option(456))

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, Seq.empty, Seq(call), new LocalBackend)
        _ <- dataAccess.updateWorkflowState(workflowId, WorkflowRunning)
        _ <- dataAccess.updateExecutionBackendInfo(workflowId, call, backendInfo)
        _ <- dataAccess.getExecutionBackendInfo(workflowId, call) map { insertResultCall =>
          insertResultCall should be(a[LocalCallBackendInfo])
          val insertResultLocalCall = insertResultCall.asInstanceOf[LocalCallBackendInfo]
          insertResultLocalCall.status should be(ExecutionStatus.Running)
          insertResultLocalCall.processId shouldNot be(empty)
          insertResultLocalCall.processId.get should be(123)
          insertResultLocalCall.resultCode shouldNot be(empty)
          insertResultLocalCall.resultCode.get should be(456)
        }
      } yield ()).futureValue
    }
  }

}
