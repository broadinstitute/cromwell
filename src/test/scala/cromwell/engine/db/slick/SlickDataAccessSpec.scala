package cromwell.engine.db.slick

import java.sql.SQLException
import java.util.UUID

import cromwell.CromwellTestkitSpec
import cromwell.CromwellTestkitSpec.TestWorkflowManagerSystem
import cromwell.binding._
import cromwell.binding.expression.WdlStandardLibraryFunctions
import cromwell.binding.types.{WdlType, WdlArrayType, WdlStringType}
import cromwell.binding.values.{WdlArray, WdlString}
import cromwell.engine.ExecutionIndex.ExecutionIndex
import cromwell.engine._
import cromwell.engine.backend.local.{LocalBackend, LocalBackendCall}
import cromwell.engine.backend.{Backend, CallLogs}
import cromwell.engine.db.{CallStatus, ExecutionDatabaseKey, LocalCallBackendInfo}
import cromwell.engine.workflow.{CallKey, WorkflowOptions}
import cromwell.parser.BackendType
import cromwell.util.SampleWdl
import cromwell.webservice.{WorkflowQueryKey, WorkflowQueryParameters}
import org.scalactic.StringNormalizations._
import org.scalatest.concurrent.ScalaFutures
import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.time.{Millis, Seconds, Span}
import org.scalatest.{FlatSpec, Matchers}

import scala.concurrent.{ExecutionContext, Future}

class SlickDataAccessSpec extends FlatSpec with Matchers with ScalaFutures {

  import TableDrivenPropertyChecks._

  val workflowManagerSystem = new TestWorkflowManagerSystem

  implicit val ec = ExecutionContext.global

  implicit val defaultPatience = PatienceConfig(timeout = Span(5, Seconds), interval = Span(100, Millis))

  lazy val localBackend = LocalBackend(workflowManagerSystem.actorSystem)
  
  val testSources = WorkflowSourceFiles("workflow test {}", "{}", "{}")

  val test2Sources = WorkflowSourceFiles("workflow test2 {}", "{}", "{}")

  object UnknownBackend extends Backend {
    type BackendCall = LocalBackendCall

    override val actorSystem = workflowManagerSystem.actorSystem

    override def adjustInputPaths(callKey: CallKey, inputs: CallInputs, workflowDescriptor: WorkflowDescriptor) =
      throw new NotImplementedError

    override def adjustOutputPaths(call: Call, outputs: CallOutputs): CallOutputs =
      throw new NotImplementedError

    override def stdoutStderr(descriptor: WorkflowDescriptor, callName: String, index: ExecutionIndex): CallLogs =
      throw new NotImplementedError

    override def initializeForWorkflow(workflow: WorkflowDescriptor) =
      throw new NotImplementedError

    override def prepareForRestart(restartableWorkflow: WorkflowDescriptor)(implicit ec: ExecutionContext) =
      throw new NotImplementedError

    override def bindCall(workflowDescriptor: WorkflowDescriptor,
                          key: CallKey,
                          locallyQualifiedInputs: CallInputs,
                          abortRegistrationFunction: AbortRegistrationFunction): BackendCall =
      throw new NotImplementedError

    override def engineFunctions(interface: IOInterface): WdlStandardLibraryFunctions =
      throw new NotImplementedError

    override def backendType: BackendType =
      throw new NotImplementedError

    override def ioInterface(workflowOptions: WorkflowOptions): IOInterface = throw new NotImplementedError
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

    /**
     *  Assert that the specified workflow has the appropriate number of calls in appropriately terminal states per
     *  the `ExpectedTerminal` function.  This function maps from an FQN + index pair to a Boolean expectation.
     */
    type ExpectedTerminal = (FullyQualifiedName, Int) => Boolean
    def assertTerminalExecution(id: WorkflowId, expectedCount: Int, expectedTerminal: ExpectedTerminal): Future[Unit] = {
      for {
        _ <- dataAccess.getExecutions(id) map { executions =>
          executions should have size expectedCount
          executions foreach { execution =>
            execution.endDt.isDefined shouldBe expectedTerminal(execution.callFqn, execution.index)
          }
        }
      } yield ()
    }

    it should "(if hsqldb) have transaction isolation mvcc" in {
      assume(canConnect || testRequired)
      import dataAccess.dataAccess.driver.api._

      val getProduct = SimpleDBIO[String](_.connection.getMetaData.getDatabaseProductName)
      //noinspection SqlDialectInspection
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
      val workflowInfo = WorkflowDescriptor(workflowId, testSources)

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, Nil, Nil, localBackend)
        _ <- dataAccess.getWorkflowsByState(Seq(WorkflowSubmitted)) map { results =>
          results shouldNot be(empty)

          val workflowResultOption = results.find(_.id == workflowId)
          workflowResultOption shouldNot be(empty)
          val workflowResult = workflowResultOption.get
          workflowResult.id should be(workflowId)
          workflowResult.sourceFiles.wdlSource should be("workflow test {}")
          workflowResult.sourceFiles.inputsJson should be("{}")
        }
      } yield ()).futureValue
    }

    it should "create and query a workflow" in {
      assume(canConnect || testRequired)
      val workflowInfo = WorkflowDescriptor(WorkflowId(UUID.randomUUID()), testSources)
      val workflow2Info = WorkflowDescriptor(WorkflowId(UUID.randomUUID()), test2Sources)

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, Nil, Nil, localBackend)
        // Update to a terminal state
        _ <- dataAccess.updateWorkflowState(workflowInfo.id, WorkflowSucceeded)
        // Put a bit of space between the two workflows
        _ = Thread.sleep(50)
        _ <- dataAccess.createWorkflow(workflow2Info, Nil, Nil, localBackend)
        // Query with no filters
        (test, test2) <- dataAccess.queryWorkflows(WorkflowQueryParameters(Seq.empty)) map { response =>
          val test = response.results find { r => r.name == "test" && r.end.isDefined } getOrElse fail
          val test2 = response.results find { _.name == "test2" } getOrElse fail
          (test, test2)
        }
        // Filter by name
        _ <- dataAccess.queryWorkflows(WorkflowQueryParameters(Seq(WorkflowQueryKey.Name.name -> "test"))) map { response =>
          val resultsByName = response.results groupBy { _.name }
          resultsByName.keys.toSet should equal(Set("test"))
        }
        // Filter by multiple names
        _ <- dataAccess.queryWorkflows(WorkflowQueryParameters(Seq(WorkflowQueryKey.Name.name -> "test", WorkflowQueryKey.Name.name -> "test2"))) map { response =>
          val resultsByName = response.results groupBy { _.name }
          resultsByName.keys.toSet should equal(Set("test", "test2"))
        }
        // Filter by status
        _ <- dataAccess.queryWorkflows(WorkflowQueryParameters(Seq(WorkflowQueryKey.Status.name -> "Submitted"))) map { response =>
          val resultsByStatus = response.results groupBy(_.status)
          resultsByStatus.keys.toSet should equal(Set("Submitted"))
        }
        // Filter by multiple statuses
        _ <- dataAccess.queryWorkflows(WorkflowQueryParameters(Seq(WorkflowQueryKey.Status.name -> "Submitted", WorkflowQueryKey.Status.name -> "Succeeded"))) map { response =>
          val resultsByStatus = response.results groupBy(_.status)
          resultsByStatus.keys.toSet should equal(Set("Submitted", "Succeeded"))
        }
        // Filter by start date
        _ <- dataAccess.queryWorkflows(WorkflowQueryParameters(Seq(WorkflowQueryKey.StartDate.name -> test2.start.toString))) map { response =>
          response.results partition { _.start.compareTo(test2.start) >= 0 } match {
            case (y, n) if y.nonEmpty && n.isEmpty => // good
            case (y, n) => fail(s"Found ${y.size} later workflows and ${n.size} earlier")
          }
        }
        // Filter by end date
        _ <- dataAccess.queryWorkflows(WorkflowQueryParameters(Seq(WorkflowQueryKey.EndDate.name -> test.end.get.toString))) map { response =>
          response.results partition { r => r.end.isDefined && r.end.get.compareTo(test.end.get) <= 0 } match {
            case (y, n) if y.nonEmpty && n.isEmpty => // good
            case (y, n) => fail(s"Found ${y.size} earlier workflows and ${n.size} later")
          }
        }
      } yield ()).futureValue
    }

    it should "query a single execution status" in {
      assume(canConnect || testRequired)
      val workflowId = WorkflowId(UUID.randomUUID())
      val workflowInfo = WorkflowDescriptor(workflowId, testSources)
      val task = new Task("taskName", Nil, Nil, Nil, null, BackendType.LOCAL)
      val callFqn = "fully.qualified.name"
      val call = new Call(None, callFqn, task, Set.empty[FullyQualifiedName], Map.empty, None)

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, Nil, Seq(call), localBackend)
        _ <- dataAccess.updateWorkflowState(workflowId, WorkflowRunning)
        _ <- dataAccess.getExecutionStatus(workflowId, ExecutionDatabaseKey(callFqn, None)) map { status =>
          status.get.executionStatus shouldBe ExecutionStatus.NotStarted
          status.get.returnCode shouldBe None
        }
      } yield ()).futureValue
    }

    it should "create and retrieve 3step.wdl with a 10,000 char pattern" in {
      assume(canConnect || testRequired)
      val workflowId = WorkflowId(UUID.randomUUID())
      val sampleWdl = SampleWdl.ThreeStepLargeJson
      val workflowInfo = WorkflowDescriptor(workflowId, sampleWdl.asWorkflowSources())

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, Nil, Nil, localBackend)
        _ <- dataAccess.getWorkflowsByState(Seq(WorkflowSubmitted)) map { results =>
          results shouldNot be(empty)

          val workflowResultOption = results.find(_.id == workflowId)
          workflowResultOption shouldNot be(empty)
          val workflowResult = workflowResultOption.get
          workflowResult.id should be(workflowId)
          workflowResult.sourceFiles.wdlSource should be(sampleWdl.wdlSource())
          workflowResult.sourceFiles.inputsJson should be(sampleWdl.wdlJson)
        }
      } yield ()).futureValue
    }

    it should "fail when saving a workflow twice" in {
      assume(canConnect || testRequired)
      val workflowId = WorkflowId(UUID.randomUUID())
      val workflowInfo = WorkflowDescriptor(workflowId, testSources)

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, Nil, Nil, localBackend)
        _ <- dataAccess.createWorkflow(workflowInfo, Nil, Nil, localBackend)
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
      val workflowInfo = WorkflowDescriptor(workflowId, testSources)

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, Nil, Nil, localBackend)
        _ <- dataAccess.updateWorkflowState(workflowId, WorkflowRunning)
        _ <- dataAccess.getWorkflowsByState(Seq(WorkflowRunning)) map { results =>
          results shouldNot be(empty)

          val workflowResultOption = results.find(_.id == workflowId)
          workflowResultOption shouldNot be(empty)
          val workflowResult = workflowResultOption.get
          workflowResult.id should be(workflowId)
          workflowResult.sourceFiles.wdlSource should be("workflow test {}")
          workflowResult.sourceFiles.inputsJson should be("{}")
        }
      } yield ()).futureValue
    }

    it should "get workflow state" in {
      assume(canConnect || testRequired)
      val workflowId = WorkflowId(UUID.randomUUID())
      val workflowInfo = WorkflowDescriptor(workflowId, testSources)
      val key = new SymbolStoreKey("myScope", "myName", None, input = true)
      val entry = new SymbolStoreEntry(key, WdlStringType, Option(new WdlString("testStringValue")))

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, Seq(entry), Nil, localBackend)
        _ <- dataAccess.getWorkflowExecution(workflowId) map { w =>
          w.status shouldBe WorkflowSubmitted.toString
          w.endDt shouldBe None
        }
        _ <- dataAccess.updateWorkflowState(workflowId, WorkflowSucceeded)
        _ <- dataAccess.getWorkflowExecution(workflowId) map { w =>
          w.status shouldBe WorkflowSucceeded.toString
          w.endDt should not be None
        }
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
        val workflowInfo = WorkflowDescriptor(workflowId, testSources)

        val task = new Task("taskName", Nil, Nil, Nil, null, BackendType.LOCAL)
        val call =
          if (setCallParent) {
            val workflow = new Workflow("workflow.name", Nil, Nil)
            val call = new Call(callAlias, "call.name", task, Set.empty[FullyQualifiedName], Map.empty, Option(workflow))
            workflow.children = Seq(call)
            call
          } else {
            new Call(callAlias, "call.name", task, Set.empty[FullyQualifiedName], Map.empty, None)
          }

        val callKey = ExecutionDatabaseKey(call.fullyQualifiedName, None)
        val shardKey = ExecutionDatabaseKey(call.fullyQualifiedName, Option(0))

        def optionallyUpdateExecutionStatus() =
          if (updateStatus)
            dataAccess.setStatus(workflowId, Seq(callKey), CallStatus(ExecutionStatus.Running, None))
          else
            Future.successful(())

        (for {
          _ <- dataAccess.createWorkflow(workflowInfo, Nil, Seq(call), localBackend)
          _ <- dataAccess.updateWorkflowState(workflowId, WorkflowRunning)
          _ <- optionallyUpdateExecutionStatus()
          _ <- dataAccess.getExecutionStatuses(workflowId) map { result =>
            result.size should be(1)
            val (key, status) = result.head
            key.fqn should be(expectedFqn)
            key.index should be(None)
            status.executionStatus should be(if (updateStatus) ExecutionStatus.Running else ExecutionStatus.NotStarted)
            status.returnCode should be(None)
          }
          _ <- dataAccess.insertCalls(workflowId, Seq(CallKey(call, Option(0))), localBackend)
          _ <- dataAccess.setStatus(workflowId, Seq(shardKey), CallStatus(ExecutionStatus.Done, Option(0)))
          _ <- dataAccess.getExecutionStatuses(workflowId) map { result =>
            result.size should be(2)
            //Previous call status should not have changed
            result should contain key callKey
            val callStatus = result.get(callKey).get
            callStatus.executionStatus should be(if (updateStatus) ExecutionStatus.Running else ExecutionStatus.NotStarted)
            callStatus.returnCode should be(None)

            result should contain key shardKey
            val shardStatus = result.get(shardKey).get
            shardStatus.executionStatus should be(ExecutionStatus.Done)
            shardStatus.returnCode should be(Option(0))
          }
          // Two calls, the call with an index is supposed to be Done and should therefore get an end date.
          _ <- assertTerminalExecution(workflowId, expectedCount = 2, expectedTerminal = (_, index) => index == 0)
        } yield ()).futureValue
      }
    }

    it should "insert a call in execution table" in {
      assume(canConnect || testRequired)
      val callFqn = "call.fully.qualified.scope"
      val workflowId = WorkflowId(UUID.randomUUID())
      val workflowInfo = WorkflowDescriptor(workflowId, testSources)
      val task = new Task("taskName", Nil, Nil, Nil, null, BackendType.LOCAL)
      val call = new Call(None, callFqn, task, Set.empty[FullyQualifiedName], Map.empty, None)
      val callKey = CallKey(call, None)
      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, Nil, Nil, localBackend)
        _ <- dataAccess.updateWorkflowState(workflowId, WorkflowRunning)
        _ <- dataAccess.insertCalls(workflowId, Seq(callKey), localBackend)
        _ <- dataAccess.getExecutionStatuses(workflowId) map { result =>
          result.size should be(1)
          val (key, status) = result.head
          key.fqn should be("call.fully.qualified.scope")
          key.index should be(None)
          status.executionStatus should be(ExecutionStatus.NotStarted)
          status.returnCode shouldBe None
        }
        _ <- assertTerminalExecution(workflowId, expectedCount = 1, expectedTerminal = (_, _) => false)
      } yield ()).futureValue
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
      val workflowInfo = WorkflowDescriptor(workflowId, testSources)
      val key = new SymbolStoreKey(callFqn, symbolFqn, None, input = true)
      val entry = new SymbolStoreEntry(key, WdlStringType, Option(new WdlString("testStringValue")))
      val task = new Task("taskName", Nil, Nil, Nil, null, BackendType.LOCAL)
      val call = new Call(None, callFqn, task, Set.empty[FullyQualifiedName], Map.empty, None)

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, Seq(entry), Nil, localBackend)
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
      val workflowInfo = WorkflowDescriptor(workflowId, testSources)
      val key = new SymbolStoreKey(callFqn, symbolFqn, None, input = true)
      val entry = new SymbolStoreEntry(key, WdlArrayType(WdlStringType), Option(wdlArray))
      val task = new Task("taskName", Nil, Nil, Nil, null, BackendType.LOCAL)
      val call = new Call(None, callFqn, task, Set.empty[FullyQualifiedName], Map.empty, None)

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, Seq(entry), Nil, localBackend)
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
      val workflowId: WorkflowId = WorkflowId(UUID.randomUUID())
      val workflowInfo = WorkflowDescriptor(workflowId, testSources)

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, Nil, Nil, localBackend)
        _ <- dataAccess.updateWorkflowState(workflowId, WorkflowRunning)
        _ <- dataAccess.getInputs(workflowId, null)
      } yield ()).failed.futureValue should be(an[IllegalArgumentException])
    }

    it should "set and get an output" in {
      assume(canConnect || testRequired)
      val callFqn = "call.fully.qualified.scope"
      val symbolLqn = "symbol"
      val workflowId = WorkflowId(UUID.randomUUID())
      val workflowInfo = WorkflowDescriptor(workflowId, testSources)
      val task = new Task("taskName", Nil, Nil, Nil, null, BackendType.LOCAL)
      val call = new Call(None, callFqn, task, Set.empty[FullyQualifiedName], Map.empty, None)
      val workflowOutputList = Seq(ReportableSymbol("call.fully.qualified.scope.symbol", WdlArrayType(WdlStringType)))

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, Nil, Nil, localBackend)
        _ <- dataAccess.updateWorkflowState(workflowId, WorkflowRunning)
        _ <- dataAccess.setOutputs(workflowId, CallKey(call, None), Map(symbolLqn -> new WdlString("testStringValue")), workflowOutputList)
        _ <- dataAccess.setOutputs(workflowId, CallKey(call, Option(0)), Map(symbolLqn -> new WdlString("testStringValueShard")), workflowOutputList)
        _ <- dataAccess.getWorkflowOutputs(workflowId) map { results =>
          results.size should be(1) //getOutputs on a workflowId does NOT return shards outputs

          val output = results.head
          output.key.scope should be("call.fully.qualified.scope")
          output.key.name should be("symbol")
          output.key.input should be(right = false) // IntelliJ highlighting
          output.wdlType should be(WdlStringType)
          output.wdlValue shouldNot be(empty)
          output.key.index should be(None)
          output.wdlValue.get should be(new WdlString("testStringValue"))
        }
      } yield ()).futureValue
    }

    it should "set and get shard statuses" in {
      assume(canConnect || testRequired)
      val callFqn1 = "call.fully.qualified.scope$s1"
      val callFqn2 = "call.fully.qualified.scope$s2"
      val workflowId = WorkflowId(UUID.randomUUID())
      val workflowInfo = WorkflowDescriptor(workflowId, testSources)
      val task = new Task("taskName", Nil, Nil, Nil, null, BackendType.LOCAL)
      val call1 = new Call(None, callFqn1, task, Set.empty[FullyQualifiedName], Map.empty, None)
      val shardIndex1 = Option(0)
      val pid1 = Option(123)
      val backendInfo1 = new LocalCallBackendInfo(pid1)
      val call2 = new Call(None, callFqn2, task, Set.empty[FullyQualifiedName], Map.empty, None)
      val shardIndex2 = Option(1)
      val pid2 = Option(987)
      val backendInfo2 = new LocalCallBackendInfo(pid2)

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, Nil, Nil, localBackend)
        _ <- dataAccess.insertCalls(workflowId, Seq(CallKey(call1, shardIndex1)), localBackend)
        _ <- dataAccess.insertCalls(workflowId, Seq(CallKey(call2, shardIndex2)), localBackend)
        _ <- dataAccess.updateExecutionBackendInfo(workflowId, CallKey(call1, shardIndex1), backendInfo1)
        _ <- dataAccess.updateExecutionBackendInfo(workflowId, CallKey(call2, shardIndex2), backendInfo2)
        _ <- dataAccess.getExecutionBackendInfo(workflowId, call1) map {
          case LocalCallBackendInfo(processId: Option[Int]) => assertResult(pid1) { processId }
          case _ => fail("Unexpected backend info type returned")
        }
        _ <- dataAccess.getExecutionBackendInfo(workflowId, call2) map {
          case LocalCallBackendInfo(processId: Option[Int]) => assertResult(pid2) { processId }
          case _ => fail("Unexpected backend info type returned")
        }
      } yield ()).futureValue
    }

    it should "set and get an output by call" in {
      assume(canConnect || testRequired)
      val callFqn = "call.fully.qualified.scope"
      val symbolLqn = "symbol"
      val workflowId = WorkflowId(UUID.randomUUID())
      val workflowInfo = WorkflowDescriptor(workflowId, testSources)
      val task = new Task("taskName", Nil, Nil, Nil, null, BackendType.LOCAL)
      val call = new Call(None, callFqn, task, Set.empty[FullyQualifiedName], Map.empty, None)

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, Nil, Nil, localBackend)
        _ <- dataAccess.updateWorkflowState(workflowId, WorkflowRunning)
        _ <- dataAccess.setOutputs(workflowId, CallKey(call, None), Map(symbolLqn -> new WdlString("testStringValue")), Seq.empty)
        _ <- dataAccess.setOutputs(workflowId, CallKey(call, Option(0)), Map(symbolLqn -> new WdlString("testStringValueShard")), Seq.empty)
        callOutput <- dataAccess.getOutputs(workflowId, ExecutionDatabaseKey(call.fullyQualifiedName, None)) map { results =>
          results.head.key.index should be(None)
          results.head.wdlValue.get should be(new WdlString("testStringValue"))
          results
        }
        shardOutput <- dataAccess.getOutputs(workflowId, ExecutionDatabaseKey(call.fullyQualifiedName, Option(0))) map { results =>
          results.head.key.index should be(Some(0))
          results.head.wdlValue.get should be(new WdlString("testStringValueShard"))
          results
        }
        _ <- Future(Seq(callOutput, shardOutput) map { results =>
          results.size should be(1)
          val resultSymbol = results.head
          val resultSymbolStoreKey = resultSymbol.key
          resultSymbolStoreKey.scope should be("call.fully.qualified.scope")
          resultSymbolStoreKey.name should be("symbol")
          resultSymbolStoreKey.input should be(right = false) // IntelliJ highlighting
          resultSymbol.wdlType should be(WdlStringType)
          resultSymbol.wdlValue shouldNot be(empty)
          results
        })
      } yield ()).futureValue
    }

    it should "fail to get outputs for a null call" in {
      assume(canConnect || testRequired)
      val workflowId = WorkflowId(UUID.randomUUID())
      val workflowInfo = WorkflowDescriptor(workflowId, testSources)

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, Nil, Nil, localBackend)
        _ <- dataAccess.updateWorkflowState(workflowId, WorkflowRunning)
        _ <- dataAccess.getOutputs(workflowId, null)
      } yield ()).failed.futureValue should be(an[IllegalArgumentException])
    }

    it should "fail to create workflow for an unknown backend" in {
      assume(canConnect || testRequired)
      val callFqn = "call.fully.qualified.scope"
      val workflowId = WorkflowId(UUID.randomUUID())
      val workflowInfo = WorkflowDescriptor(workflowId, testSources)
      val task = new Task("taskName", Nil, Nil, Nil, null, BackendType.LOCAL)
      val call = new Call(None, callFqn, task, Set.empty[FullyQualifiedName], Map.empty, None)

      dataAccess.createWorkflow(workflowInfo, Nil, Seq(call),
        UnknownBackend).failed.futureValue should be(an[IllegalArgumentException])
    }

    it should "fail to create workflow for a null backend" in {
      assume(canConnect || testRequired)
      val callFqn = "call.fully.qualified.scope"
      val workflowId = WorkflowId(UUID.randomUUID())
      val workflowInfo = WorkflowDescriptor(workflowId, testSources)
      val task = new Task("taskName", Nil, Nil, Nil, null, BackendType.LOCAL)
      val call = new Call(None, callFqn, task, Set.empty[FullyQualifiedName], Map.empty, None)

      dataAccess.createWorkflow(workflowInfo, Nil, Seq(call),
        null).failed.futureValue should be(an[IllegalArgumentException])
    }

    it should "set and get the same symbol with IO as input then output" in {
      assume(canConnect || testRequired)
      val callFqn = "call.fully.qualified.scope"
      val symbolLqn = "symbol"
      val symbolFqn = callFqn + "." + symbolLqn
      val workflowId = WorkflowId(UUID.randomUUID())
      val workflowInfo = WorkflowDescriptor(workflowId, testSources)
      val key = new SymbolStoreKey(callFqn, symbolFqn, None, input = true)
      val entry = new SymbolStoreEntry(key, WdlStringType, Option(new WdlString("testStringValue")))
      val task = new Task("taskName", Nil, Nil, Nil, null, BackendType.LOCAL)
      val call = new Call(None, callFqn, task, Set.empty[FullyQualifiedName], Map.empty, None)

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, Seq(entry), Nil, localBackend)
        _ <- dataAccess.updateWorkflowState(workflowId, WorkflowRunning)
        _ <- dataAccess.setOutputs(workflowId, CallKey(call, None), Map(symbolLqn -> new WdlString("testStringValue")), Seq.empty)
        _ <- dataAccess.getOutputs(workflowId, ExecutionDatabaseKey(call.fullyQualifiedName, None)) map { results =>
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
      val workflowInfo = WorkflowDescriptor(workflowId, testSources)
      val task = new Task("taskName", Nil, Nil, Nil, null, BackendType.LOCAL)
      val call = new Call(None, callFqn, task, Set.empty[FullyQualifiedName], Map.empty, None)

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, Seq(), Nil, localBackend)
        _ <- dataAccess.updateWorkflowState(workflowId, WorkflowRunning)
        _ <- dataAccess.setOutputs(workflowId, CallKey(call, None), Map(symbolFqn -> new WdlString("testStringValue")), Seq.empty)
        // This second attempt should now fail:
        _ <- dataAccess.setOutputs(workflowId, CallKey(call, None), Map(symbolFqn -> new WdlString("testStringValue")), Seq.empty)
      } yield ()).failed.futureValue should be(a[SQLException])
    }

    it should "set and get a backend info" in {
      assume(canConnect || testRequired)
      val workflowId = WorkflowId(UUID.randomUUID())
      val workflowInfo = WorkflowDescriptor(workflowId, testSources)
      val task = new Task("taskName", Nil, Nil, Nil, null, BackendType.LOCAL)
      val call = new Call(None, "fully.qualified.name", task, Set.empty[FullyQualifiedName], Map.empty, None)
      val backendInfo = new LocalCallBackendInfo(Option(123))

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, Nil, Seq(call), localBackend)
        _ <- dataAccess.updateExecutionBackendInfo(workflowId, CallKey(call, None), backendInfo)
        _ <- dataAccess.getExecutionBackendInfo(workflowId, call) map { insertResultCall =>
          insertResultCall should be(a[LocalCallBackendInfo])
          val insertResultLocalCall = insertResultCall.asInstanceOf[LocalCallBackendInfo]
          insertResultLocalCall.processId shouldNot be(empty)
          insertResultLocalCall.processId.get should be(123)
        }
      } yield ()).futureValue
    }

    // Queries use `.head` a lot. There was a bug that pulled the backend info by fqn, but for any workflow.
    it should "set and get a backend info for same call on two workflows" in {
      assume(canConnect || testRequired)
      val workflowId1 = WorkflowId(UUID.randomUUID())
      val workflowId2 = WorkflowId(UUID.randomUUID())
      val workflowInfo1 = WorkflowDescriptor(workflowId1, testSources)
      val workflowInfo2 = WorkflowDescriptor(workflowId2, testSources)
      val task = new Task("taskName", Nil, Nil, Nil, null, BackendType.LOCAL)
      val call = new Call(None, "fully.qualified.name", task, Set.empty[FullyQualifiedName], Map.empty, None)
      val backendInfo1 = new LocalCallBackendInfo(Option(123))
      val backendInfo2 = new LocalCallBackendInfo(Option(321))

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo1, Nil, Seq(call), localBackend)
        _ <- dataAccess.createWorkflow(workflowInfo2, Nil, Seq(call), localBackend)
        _ <- dataAccess.updateExecutionBackendInfo(workflowId1, CallKey(call, None), backendInfo1)
        _ <- dataAccess.updateExecutionBackendInfo(workflowId2, CallKey(call, None), backendInfo2)
        _ <- dataAccess.getExecutionBackendInfo(workflowId1, call) map { insertResultCall =>
          insertResultCall should be(a[LocalCallBackendInfo])
          val insertResultLocalCall = insertResultCall.asInstanceOf[LocalCallBackendInfo]
          insertResultLocalCall.processId shouldNot be(empty)
          insertResultLocalCall.processId.get should be(123)
        }
        _ <- dataAccess.getExecutionBackendInfo(workflowId2, call) map { insertResultCall =>
          insertResultCall should be(a[LocalCallBackendInfo])
          val insertResultLocalCall = insertResultCall.asInstanceOf[LocalCallBackendInfo]
          insertResultLocalCall.processId shouldNot be(empty)
          insertResultLocalCall.processId.get should be(321)
        }
      } yield ()).futureValue
    }

    it should "fail to set null backend info" in {
      assume(canConnect || testRequired)
      val workflowId = WorkflowId(UUID.randomUUID())
      val workflowInfo = WorkflowDescriptor(workflowId, testSources)
      val task = new Task("taskName", Nil, Nil, Nil, null, BackendType.LOCAL)
      val call = new Call(None, "fully.qualified.name", task, Set.empty[FullyQualifiedName], Map.empty, None)

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, Nil, Seq(call), localBackend)
        _ <- dataAccess.updateWorkflowState(workflowId, WorkflowRunning)
        _ <- dataAccess.updateExecutionBackendInfo(workflowId, CallKey(call, None), null)
      } yield ()).failed.futureValue should be(an[IllegalArgumentException])
    }

    it should "update call start and end dates appropriately" in {
      assume(canConnect || testRequired)
      val callFqn = "call.fully.qualified.scope"
      val workflowId = WorkflowId(UUID.randomUUID())
      val workflowInfo = WorkflowDescriptor(workflowId, testSources)
      val task = new Task("taskName", Nil, Nil, Nil, null, BackendType.LOCAL)
      val call = new Call(None, callFqn, task, Set.empty[FullyQualifiedName], Map.empty, None)

      // Assert a singular `Execution` in `executions`, and that the dates of the `Execution` are defined
      // or not as specified by `startDefined` and `endDefined`.
      def assertDates(startDefined: Boolean, endDefined: Boolean)(executions: Traversable[Execution]): Unit = {
        executions should have size 1
        executions foreach { e =>
          e.startDt.isDefined shouldBe startDefined
          e.endDt.isDefined shouldBe endDefined
        }
      }

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, Nil, Seq(call), localBackend)
        callKeys = Seq(ExecutionDatabaseKey(callFqn, None))

        _ <- dataAccess.setStatus(workflowId, callKeys, ExecutionStatus.NotStarted)
        _ <- dataAccess.getExecutions(workflowId) map assertDates(startDefined = false, endDefined = false)

        _ <- dataAccess.setStatus(workflowId, callKeys, ExecutionStatus.Starting)
        _ <- dataAccess.getExecutions(workflowId) map assertDates(startDefined = true, endDefined = false)

        _ <- dataAccess.setStatus(workflowId, callKeys, ExecutionStatus.Done)
        _ <- dataAccess.getExecutions(workflowId) map assertDates(startDefined = true, endDefined = true)
      } yield()).futureValue
    }

    it should "shutdown the database" in {
      assume(canConnect || testRequired)
      dataAccess.shutdown().futureValue
    }
  }
}
