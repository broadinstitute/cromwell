package cromwell.engine.db.slick

import java.sql.SQLException
import java.util.UUID

import cromwell.binding._
import cromwell.binding.command.CommandPart
import cromwell.binding.types.{WdlArrayType, WdlStringType}
import cromwell.binding.values.{WdlArray, WdlString, WdlValue}
import cromwell.engine._
import cromwell.engine.backend.Backend.RestartableWorkflow
import cromwell.engine.backend.local.{LocalBackendCall, LocalBackend}
import cromwell.engine.backend.{Backend, BackendCall, StdoutStderr}
import cromwell.engine.db.DataAccess.WorkflowInfo
import cromwell.engine.db.{DataAccess, LocalCallBackendInfo}
import cromwell.parser.BackendType
import cromwell.util.SampleWdl
import org.scalactic.StringNormalizations._
import org.scalatest.concurrent.ScalaFutures
import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.time.{Millis, Seconds, Span}
import org.scalatest.{FlatSpec, Matchers}

import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Success, Try}

class SlickDataAccessSpec extends FlatSpec with Matchers with ScalaFutures {

  import TableDrivenPropertyChecks._

  implicit val ec = ExecutionContext.global

  implicit val defaultPatience = PatienceConfig(timeout = Span(5, Seconds), interval = Span(100, Millis))

  lazy val localBackend = new LocalBackend

  object UnknownBackend extends Backend {
    type BackendCall = LocalBackendCall
    override def adjustInputPaths(call: Call, inputs: CallInputs) = Map.empty
    override def adjustOutputPaths(call: Call, outputs: CallOutputs): CallOutputs = outputs
    override def stdoutStderr(workflowId: WorkflowId, workflowName: String, callName: String): StdoutStderr = ???

    override def initializeForWorkflow(workflow: WorkflowDescriptor) = Success(Map.empty)

    override def handleCallRestarts(restartableWorkflows: Seq[RestartableWorkflow],
                                    dataAccess: DataAccess)(implicit ec: ExecutionContext) = Future.successful(())

    override def bindCall(workflowDescriptor: WorkflowDescriptor,
                               call: Call,
                               locallyQualifiedInputs: CallInputs,
                               abortRegistrationFunction: AbortRegistrationFunction): BackendCall = ???
    override def execute(bc: BackendCall): Try[Map[String, WdlValue]] = Success(Map.empty)

    override def backendType: BackendType = ???
  }

  // Tests against main database used for command line
  "SlickDataAccess (main.hsqldb)" should behave like databaseWithConfig("main.hsqldb", testRequired = true)

  // Tests using liquibase, but in memory
  "SlickDataAccess (test.hsqldb)" should behave like databaseWithConfig("test.hsqldb", testRequired = true)

  // If able to connect, then also run the tests on mysql, but it's not required
  "SlickDataAccess (test.mysql)" should behave like databaseWithConfig("test.mysql", testRequired = false)

  def databaseWithConfig(path: => String, testRequired: => Boolean): Unit = {

    lazy val testDatabase = new TestSlickDatabase(path)
    lazy val canConnect = testRequired || testDatabase.isValidConnection.futureValue
    lazy val dataAccess = testDatabase.slickDataAccess

    it should "(if hsqldb) have transaction isolation mvcc" in {
      assume(canConnect || testRequired)
      import dataAccess.dataAccess.driver.api._

      val getProduct = SimpleDBIO[String](_.connection.getMetaData.getDatabaseProductName)
      val getHsqldbTx = sql"""SELECT PROPERTY_VALUE
                              FROM INFORMATION_SCHEMA.SYSTEM_PROPERTIES
                              WHERE PROPERTY_NAME = 'hsqldb.tx'""".as[String].head

      (for {
        product <- dataAccess.database.run(getProduct)
        _ <- product match {
          case "HSQL Database Engine" =>
            dataAccess.database.run(getHsqldbTx) map { hsqldbTx =>
              (hsqldbTx shouldEqual "mvcc")(after being lowerCased)
            }
          case _ => Future.successful(())
        }
      } yield ()).futureValue
    }

    it should "setup via liquibase if necessary" in {
      assume(canConnect || testRequired)
      if (testDatabase.useLiquibase)
        testDatabase.setupLiquibase()
    }

    it should "create and retrieve the workflow for just reading" in {
      assume(canConnect || testRequired)
      val workflowId = WorkflowId(UUID.randomUUID())
      val workflowInfo = new WorkflowInfo(workflowId, "source", "{}")

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, Seq.empty, Seq.empty, localBackend)
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

    it should "query a single execution status" in {
      assume(canConnect || testRequired)
      val workflowId = WorkflowId(UUID.randomUUID())
      val workflowInfo = new WorkflowInfo(workflowId, "source", "{}")
      val task = new Task("taskName", Seq.empty[Declaration], Seq.empty[CommandPart], Seq.empty, null, BackendType.LOCAL)
      val callFqn = "fully.qualified.name"
      val call = new Call(None, callFqn, task, Map.empty, None)
      val backendInfo = new LocalCallBackendInfo(ExecutionStatus.Running, Option(123), Option(456))

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, Seq.empty, Seq(call), localBackend)
        _ <- dataAccess.updateWorkflowState(workflowId, WorkflowRunning)
        _ <- dataAccess.getExecutionStatus(workflowId, callFqn) map { status =>
          status.get shouldBe ExecutionStatus.NotStarted
        }
      } yield ()).futureValue
    }

    it should "create and retrieve 3step.wdl with a 10,000 char pattern" in {
      assume(canConnect || testRequired)
      val workflowId = WorkflowId(UUID.randomUUID())
      val sampleWdl = SampleWdl.ThreeStepLargeJson
      val workflowInfo = new WorkflowInfo(workflowId, sampleWdl.wdlSource(), sampleWdl.wdlJson)

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, Seq.empty, Seq.empty, localBackend)
        _ <- dataAccess.getWorkflowsByState(Seq(WorkflowSubmitted)) map { results =>
          results shouldNot be(empty)

          val workflowResultOption = results.find(_.workflowId == workflowId)
          workflowResultOption shouldNot be(empty)
          val workflowResult = workflowResultOption.get
          workflowResult.workflowId should be(workflowId)
          workflowResult.wdlSource should be(sampleWdl.wdlSource())
          workflowResult.wdlJson should be(sampleWdl.wdlJson)
        }
      } yield ()).futureValue
    }

    it should "fail when saving a workflow twice" in {
      assume(canConnect || testRequired)
      val workflowId = WorkflowId(UUID.randomUUID())
      val workflowInfo = new WorkflowInfo(workflowId, "source", "{}")

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, Seq.empty, Seq.empty, localBackend)
        _ <- dataAccess.createWorkflow(workflowInfo, Seq.empty, Seq.empty, localBackend)
      } yield ()).failed.futureValue should be(a[SQLException])
    }

    it should "fail when updating a non-existent workflow state" in {
      assume(canConnect || testRequired)
      val workflowId = WorkflowId(UUID.randomUUID())

      (for {
        _ <- dataAccess.updateWorkflowState(workflowId, WorkflowRunning)
      } yield ()).failed.futureValue should be(an[IllegalArgumentException])
    }

    it should "update and get a workflow state" in {
      assume(canConnect || testRequired)
      val workflowId = WorkflowId(UUID.randomUUID())
      val workflowInfo = new WorkflowInfo(workflowId, "source", "{}")

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, Seq.empty, Seq.empty, localBackend)
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

    it should "get workflow state" in {
      assume(canConnect || testRequired)
      val workflowId = WorkflowId(UUID.randomUUID())
      val workflowInfo = new WorkflowInfo(workflowId, "source", "{}")
      val key = new SymbolStoreKey("myScope", "myName", None, input = true)
      val entry = new SymbolStoreEntry(key, WdlStringType, Option(new WdlString("testStringValue")))

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, Seq(entry), Seq.empty, localBackend)
        _ <- dataAccess.updateWorkflowState(workflowId, WorkflowSucceeded)
        _ <- dataAccess.getWorkflowState(workflowId) map { result =>
          result shouldNot be(empty)
          result.get should be(WorkflowSucceeded)
        }
      } yield ()).futureValue
    }

    // We were using call.taskFqn instead of call.fullyQualifiedName.
    // Make sure all permutations of aliases, nested calls, & updates are working correctly.
    val executionStatusPermutations = Table(
      ("updateStatus", "useAlias", "setCallParent", "expectedFqn"),
      (false, false, false, "call.name"),
      (false, true, false, "call.alias"),
      (true, false, false, "call.name"),
      (true, true, false, "call.alias"),
      (false, false, true, "workflow.name.call.name"),
      (false, true, true, "workflow.name.call.alias"),
      (true, false, true, "workflow.name.call.name"),
      (true, true, true, "workflow.name.call.alias"))

    forAll(executionStatusPermutations) { (updateStatus, useAlias, setCallParent, expectedFqn) =>

      val spec = "%s execution status for a call%s%s".format(
        if (updateStatus) "updated" else "initial",
        if (setCallParent) " in workflow" else "",
        if (useAlias) " with alias" else "")

      val callAlias = if (useAlias) Some("call.alias") else None

      it should s"get $spec" in {
        assume(canConnect || testRequired)
        val workflowId = WorkflowId(UUID.randomUUID())
        val workflowInfo = new WorkflowInfo(workflowId, "source", "{}")

        val task = new Task("taskName", Seq.empty[Declaration], Seq.empty[CommandPart], Seq.empty, null, BackendType.LOCAL)
        val call = if (setCallParent) {
          val w = new Workflow("workflow.name", Seq.empty[Declaration], None, Seq.empty[WorkflowOutputDeclaration])
          val c = new Call(callAlias, "call.name", task, Map.empty, Option(w))
          w.setChildren(Seq(c))
          c
        } else new Call(callAlias, "call.name", task, Map.empty, None)

        def optionallyUpdateExecutionStatus() =
          if (updateStatus)
            dataAccess.setStatus(workflowId, Seq(call.fullyQualifiedName), ExecutionStatus.Running)
          else
            Future.successful(())

        (for {
          _ <- dataAccess.createWorkflow(workflowInfo, Seq.empty, Seq(call), localBackend)
          _ <- dataAccess.updateWorkflowState(workflowId, WorkflowRunning)
          _ <- optionallyUpdateExecutionStatus()
          _ <- dataAccess.getExecutionStatuses(workflowId) map { result =>
            result.size should be(1)
            val (fqn, status) = result.head
            fqn should be(expectedFqn)
            status should be(if (updateStatus) ExecutionStatus.Running else ExecutionStatus.NotStarted)
          }
        } yield ()).futureValue
      }

    }

    it should "fail to get an non-existent execution status" in {
      assume(canConnect || testRequired)
      dataAccess.getExecutionStatuses(WorkflowId(UUID.randomUUID())).failed.futureValue should be(a[NoSuchElementException])
    }

    it should "get a symbol input" in {
      assume(canConnect || testRequired)
      val callFqn = "call.fully.qualified.scope"
      val symbolFqn = "symbol.fully.qualified.scope"
      val workflowId = WorkflowId(UUID.randomUUID())
      val workflowInfo = new WorkflowInfo(workflowId, "source", "{}")
      val key = new SymbolStoreKey(callFqn, symbolFqn, None, input = true)
      val entry = new SymbolStoreEntry(key, WdlStringType, Option(new WdlString("testStringValue")))
      val task = new Task("taskName", Seq.empty[Declaration], Seq.empty[CommandPart], Seq.empty, null, BackendType.LOCAL)
      val call = new Call(None, callFqn, task, Map.empty, None)

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, Seq(entry), Seq.empty, localBackend)
        _ <- dataAccess.updateWorkflowState(workflowId, WorkflowRunning)
        _ <- dataAccess.getInputs(workflowId, call) map { resultSymbols =>
          resultSymbols.size should be(1)
          val resultSymbol = resultSymbols.head
          val resultSymbolStoreKey = resultSymbol.key
          resultSymbolStoreKey.scope should be("call.fully.qualified.scope")
          resultSymbolStoreKey.name should be("symbol.fully.qualified.scope")
          resultSymbolStoreKey.index should be(None)
          resultSymbolStoreKey.input should be(right = true) // IntelliJ highlighting
          resultSymbol.wdlType should be(WdlStringType)
          resultSymbol.wdlValue shouldNot be(empty)
          resultSymbol.wdlValue.get should be(new WdlString("testStringValue"))
        }
      } yield ()).futureValue
    }


    it should "get a symbol input that has a very long WDL value field" in {
      assume(canConnect || testRequired)
      val wdlArray = new WdlArray(WdlArrayType(WdlStringType), Seq(WdlString("test"), WdlString("*" * 10000)))
      val callFqn = "call.fully.qualified.scope"
      val symbolFqn = "symbol.fully.qualified.scope"
      val workflowId = WorkflowId(UUID.randomUUID())
      val workflowInfo = new WorkflowInfo(workflowId, "source", "{}")
      val key = new SymbolStoreKey(callFqn, symbolFqn, None, input = true)
      val entry = new SymbolStoreEntry(key, WdlArrayType(WdlStringType), Option(wdlArray))
      val task = new Task("taskName", Seq.empty[Declaration], Seq.empty[CommandPart], Seq.empty, null, BackendType.LOCAL)
      val call = new Call(None, callFqn, task, Map.empty, None)

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, Seq(entry), Seq.empty, localBackend)
        _ <- dataAccess.updateWorkflowState(workflowId, WorkflowRunning)
        _ <- dataAccess.getInputs(workflowId, call) map { resultSymbols =>
          resultSymbols.size should be(1)
          val resultSymbol = resultSymbols.head
          val resultSymbolStoreKey = resultSymbol.key
          resultSymbolStoreKey.scope should be("call.fully.qualified.scope")
          resultSymbolStoreKey.name should be("symbol.fully.qualified.scope")
          resultSymbolStoreKey.index should be(None)
          resultSymbolStoreKey.input should be(right = true) // IntelliJ highlighting
          resultSymbol.wdlType should be(WdlArrayType(WdlStringType))
          resultSymbol.wdlValue shouldNot be(empty)
          resultSymbol.wdlValue.get should be(wdlArray)
        }
      } yield ()).futureValue
    }

    it should "fail to get inputs for a null call" in {
      assume(canConnect || testRequired)
      val workflowId:WorkflowId = WorkflowId(UUID.randomUUID())
      val workflowInfo = new WorkflowInfo(workflowId, "source", "{}")

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, Seq.empty, Seq.empty, localBackend)
        _ <- dataAccess.updateWorkflowState(workflowId, WorkflowRunning)
        _ <- dataAccess.getInputs(workflowId, null)
      } yield ()).failed.futureValue should be(an[IllegalArgumentException])
    }

    it should "set and get an output" in {
      assume(canConnect || testRequired)
      val callFqn = "call.fully.qualified.scope"
      val symbolLqn = "symbol"
      val workflowId = WorkflowId(UUID.randomUUID())
      val workflowInfo = new WorkflowInfo(workflowId, "source", "{}")
      val task = new Task("taskName", Seq.empty[Declaration], Seq.empty[CommandPart], Seq.empty, null, BackendType.LOCAL)
      val call = new Call(None, callFqn, task, Map.empty, None)

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, Seq.empty, Seq.empty, localBackend)
        _ <- dataAccess.updateWorkflowState(workflowId, WorkflowRunning)
        _ <- dataAccess.setOutputs(workflowId, call, Map(symbolLqn -> new WdlString("testStringValue")))
        _ <- dataAccess.getOutputs(workflowId) map { results =>
          results.size should be(1)
          val resultSymbol = results.head
          val resultSymbolStoreKey = resultSymbol.key
          resultSymbolStoreKey.scope should be("call.fully.qualified.scope")
          resultSymbolStoreKey.name should be("symbol")
          resultSymbolStoreKey.index should be(None)
          resultSymbolStoreKey.input should be(right = false) // IntelliJ highlighting
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
      val workflowId = WorkflowId(UUID.randomUUID())
      val workflowInfo = new WorkflowInfo(workflowId, "source", "{}")
      val task = new Task("taskName", Seq.empty[Declaration], Seq.empty[CommandPart], Seq.empty, null, BackendType.LOCAL)
      val call = new Call(None, callFqn, task, Map.empty, None)

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, Seq.empty, Seq.empty, localBackend)
        _ <- dataAccess.updateWorkflowState(workflowId, WorkflowRunning)
        _ <- dataAccess.setOutputs(workflowId, call, Map(symbolLqn -> new WdlString("testStringValue")))
        _ <- dataAccess.getOutputs(workflowId, call.fullyQualifiedName) map { results =>
          results.size should be(1)
          val resultSymbol = results.head
          val resultSymbolStoreKey = resultSymbol.key
          resultSymbolStoreKey.scope should be("call.fully.qualified.scope")
          resultSymbolStoreKey.name should be("symbol")
          resultSymbolStoreKey.index should be(None)
          resultSymbolStoreKey.input should be(right = false) // IntelliJ highlighting
          resultSymbol.wdlType should be(WdlStringType)
          resultSymbol.wdlValue shouldNot be(empty)
          resultSymbol.wdlValue.get should be(new WdlString("testStringValue"))
        }
      } yield ()).futureValue
    }

    it should "fail to get outputs for a null call" in {
      assume(canConnect || testRequired)
      val workflowId = WorkflowId(UUID.randomUUID())
      val workflowInfo = new WorkflowInfo(workflowId, "source", "{}")

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, Seq.empty, Seq.empty, localBackend)
        _ <- dataAccess.updateWorkflowState(workflowId, WorkflowRunning)
        _ <- dataAccess.getOutputs(workflowId, null)
      } yield ()).failed.futureValue should be(an[IllegalArgumentException])
    }

    it should "fail to create workflow for an unknown backend" in {
      assume(canConnect || testRequired)
      val callFqn = "call.fully.qualified.scope"
      val workflowId = WorkflowId(UUID.randomUUID())
      val workflowInfo = new WorkflowInfo(workflowId, "source", "{}")
      val task = new Task("taskName", Seq.empty[Declaration], Seq.empty[CommandPart], Seq.empty, null, BackendType.LOCAL)
      val call = new Call(None, callFqn, task, Map.empty, None)

      dataAccess.createWorkflow(workflowInfo, Seq.empty, Seq(call),
        UnknownBackend).failed.futureValue should be(an[IllegalArgumentException])
    }

    it should "fail to create workflow for a null backend" in {
      assume(canConnect || testRequired)
      val callFqn = "call.fully.qualified.scope"
      val workflowId = WorkflowId(UUID.randomUUID())
      val workflowInfo = new WorkflowInfo(workflowId, "source", "{}")
      val task = new Task("taskName", Seq.empty[Declaration], Seq.empty[CommandPart], Seq.empty, null, BackendType.LOCAL)
      val call = new Call(None, callFqn, task, Map.empty, None)

      dataAccess.createWorkflow(workflowInfo, Seq.empty, Seq(call),
        null).failed.futureValue should be(an[IllegalArgumentException])
    }

    it should "set and get the same symbol with IO as input then output" in {
      assume(canConnect || testRequired)
      val callFqn = "call.fully.qualified.scope"
      val symbolLqn = "symbol"
      val symbolFqn = callFqn + "." + symbolLqn
      val workflowId = WorkflowId(UUID.randomUUID())
      val workflowInfo = new WorkflowInfo(workflowId, "source", "{}")
      val key = new SymbolStoreKey(callFqn, symbolFqn, None, input = true)
      val entry = new SymbolStoreEntry(key, WdlStringType, Option(new WdlString("testStringValue")))
      val task = new Task("taskName", Seq.empty[Declaration], Seq.empty[CommandPart], Seq.empty, null, BackendType.LOCAL)
      val call = new Call(None, callFqn, task, Map.empty, None)

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, Seq(entry), Seq.empty, localBackend)
        _ <- dataAccess.updateWorkflowState(workflowId, WorkflowRunning)
        _ <- dataAccess.setOutputs(workflowId, call, Map(symbolLqn -> new WdlString("testStringValue")))
        _ <- dataAccess.getOutputs(workflowId, call.fullyQualifiedName) map { results =>
          results.size should be(1)
          val resultSymbol = results.head
          val resultSymbolStoreKey = resultSymbol.key
          resultSymbolStoreKey.scope should be("call.fully.qualified.scope")
          resultSymbolStoreKey.name should be("symbol")
          resultSymbolStoreKey.index should be(None)
          resultSymbolStoreKey.input should be(right = false) // IntelliJ highlighting
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
      val workflowId = WorkflowId(UUID.randomUUID())
      val workflowInfo = new WorkflowInfo(workflowId, "source", "{}")
      val key = new SymbolStoreKey(callFqn, symbolFqn, None, input = false)
      val entry = new SymbolStoreEntry(key, WdlStringType, Option(new WdlString("testStringValue")))
      val task = new Task("taskName", Seq.empty[Declaration], Seq.empty[CommandPart], Seq.empty, null, BackendType.LOCAL)
      val call = new Call(None, callFqn, task, Map.empty, None)

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, Seq(entry), Seq.empty, localBackend)
        _ <- dataAccess.updateWorkflowState(workflowId, WorkflowRunning)
        _ <- dataAccess.setOutputs(workflowId, call, Map(symbolFqn -> new WdlString("testStringValue")))
      } yield ()).failed.futureValue should be(a[SQLException])
    }

    it should "set and get a backend info" in {
      assume(canConnect || testRequired)
      val workflowId = WorkflowId(UUID.randomUUID())
      val workflowInfo = new WorkflowInfo(workflowId, "source", "{}")
      val task = new Task("taskName", Seq.empty[Declaration], Seq.empty[CommandPart], Seq.empty, null, BackendType.LOCAL)
      val call = new Call(None, "fully.qualified.name", task, Map.empty, None)
      val backendInfo = new LocalCallBackendInfo(ExecutionStatus.Running, Option(123), Option(456))

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, Seq.empty, Seq(call), localBackend)
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

    // Queries use `.head` a lot. There was a bug that pulled the backend info by fqn, but for any workflow.
    it should "set and get a backend info for same call on two workflows" in {
      assume(canConnect || testRequired)
      val workflowId1 = WorkflowId(UUID.randomUUID())
      val workflowId2 = WorkflowId(UUID.randomUUID())
      val workflowInfo1 = new WorkflowInfo(workflowId1, "source", "{}")
      val workflowInfo2 = new WorkflowInfo(workflowId2, "source", "{}")
      val task = new Task("taskName", Seq.empty[Declaration], Seq.empty[CommandPart], Seq.empty, null, BackendType.LOCAL)
      val call = new Call(None, "fully.qualified.name", task, Map.empty, None)
      val backendInfo1 = new LocalCallBackendInfo(ExecutionStatus.Running, Option(123), Option(456))
      val backendInfo2 = new LocalCallBackendInfo(ExecutionStatus.Failed, Option(321), Option(654))

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo1, Seq.empty, Seq(call), localBackend)
        _ <- dataAccess.createWorkflow(workflowInfo2, Seq.empty, Seq(call), localBackend)
        _ <- dataAccess.updateWorkflowState(workflowId1, WorkflowRunning)
        _ <- dataAccess.updateWorkflowState(workflowId2, WorkflowRunning)
        _ <- dataAccess.updateExecutionBackendInfo(workflowId1, call, backendInfo1)
        _ <- dataAccess.updateExecutionBackendInfo(workflowId2, call, backendInfo2)
        _ <- dataAccess.getExecutionBackendInfo(workflowId1, call) map { insertResultCall =>
          insertResultCall should be(a[LocalCallBackendInfo])
          val insertResultLocalCall = insertResultCall.asInstanceOf[LocalCallBackendInfo]
          insertResultLocalCall.status should be(ExecutionStatus.Running)
          insertResultLocalCall.processId shouldNot be(empty)
          insertResultLocalCall.processId.get should be(123)
          insertResultLocalCall.resultCode shouldNot be(empty)
          insertResultLocalCall.resultCode.get should be(456)
        }
        _ <- dataAccess.getExecutionBackendInfo(workflowId2, call) map { insertResultCall =>
          insertResultCall should be(a[LocalCallBackendInfo])
          val insertResultLocalCall = insertResultCall.asInstanceOf[LocalCallBackendInfo]
          insertResultLocalCall.status should be(ExecutionStatus.Failed)
          insertResultLocalCall.processId shouldNot be(empty)
          insertResultLocalCall.processId.get should be(321)
          insertResultLocalCall.resultCode shouldNot be(empty)
          insertResultLocalCall.resultCode.get should be(654)
        }
      } yield ()).futureValue
    }

    it should "fail to set null backend info" in {
      assume(canConnect || testRequired)
      val workflowId = WorkflowId(UUID.randomUUID())
      val workflowInfo = new WorkflowInfo(workflowId, "source", "{}")
      val task = new Task("taskName", Seq.empty[Declaration], Seq.empty[CommandPart], Seq.empty, null, BackendType.LOCAL)
      val call = new Call(None, "fully.qualified.name", task, Map.empty, None)

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, Seq.empty, Seq(call), localBackend)
        _ <- dataAccess.updateWorkflowState(workflowId, WorkflowRunning)
        _ <- dataAccess.updateExecutionBackendInfo(workflowId, call, null)
      } yield ()).failed.futureValue should be(an[IllegalArgumentException])
    }

    it should "shutdown the database" in {
      assume(canConnect || testRequired)
      dataAccess.shutdown().futureValue
    }
  }
}
