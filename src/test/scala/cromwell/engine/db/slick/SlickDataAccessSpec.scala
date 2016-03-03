package cromwell.engine.db.slick

import java.sql.SQLException
import java.util.UUID

import better.files._
import com.google.api.client.util.ExponentialBackOff
import com.typesafe.config.ConfigFactory
import cromwell.CromwellSpec.IntegrationTest
import cromwell.CromwellTestkitSpec.TestWorkflowManagerSystem
import cromwell.engine.Hashing._
import cromwell.engine.backend.local.LocalBackend.InfoKeys
import cromwell.engine.backend.local.{LocalBackend, LocalBackendCall}
import cromwell.engine.backend.{Backend, BackendType, CallLogs}
import cromwell.engine.db.slick.SlickDataAccessSpec.{AllowFalse, AllowTrue}
import cromwell.engine.db.{DiffResultFilter, ExecutionDatabaseKey}
import cromwell.engine.io.IoInterface
import cromwell.engine.workflow.{BackendCallKey, ScatterKey, WorkflowOptions}
import cromwell.engine.{CallOutput, CallOutputs, _}
import cromwell.util.SampleWdl
import cromwell.webservice
import cromwell.webservice.{CallCachingParameters, WorkflowQueryKey, WorkflowQueryParameters}
import org.joda.time.DateTime
import org.scalactic.StringNormalizations._
import org.scalatest.PartialFunctionValues._
import org.scalatest.concurrent.PatienceConfiguration.Timeout
import org.scalatest.concurrent.ScalaFutures
import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.time.SpanSugar._
import org.scalatest.time.{Millis, Seconds, Span}
import org.scalatest.{BeforeAndAfterAll, FlatSpec, Matchers}
import org.specs2.mock.Mockito
import wdl4s.types.{WdlArrayType, WdlStringType}
import wdl4s.values.{WdlArray, WdlString}
import wdl4s.{CallInputs, _}

import scala.concurrent.{ExecutionContext, Future}

object SlickDataAccessSpec {
  val AllowFalse = Seq(webservice.QueryParameter("allow", "false"))
  val AllowTrue = Seq(webservice.QueryParameter("allow", "true"))
}

class SlickDataAccessSpec extends FlatSpec with Matchers with ScalaFutures with BeforeAndAfterAll with Mockito {

  import TableDrivenPropertyChecks._

  val workflowManagerSystem = new TestWorkflowManagerSystem

  override protected def afterAll() = {
    workflowManagerSystem.shutdownTestActorSystem()
    super.afterAll()
  }

  implicit val ec = ExecutionContext.global

  implicit val defaultPatience = PatienceConfig(timeout = Span(5, Seconds), interval = Span(100, Millis))

  lazy val localBackend = LocalBackend(workflowManagerSystem.actorSystem)

  val testSources = WorkflowSourceFiles("workflow test {}", "{}", "{}")

  implicit val hasher = WorkflowDescriptor(WorkflowId(UUID.randomUUID()), testSources).fileHasher

  val test2Sources = WorkflowSourceFiles("workflow test2 {}", "{}", "{}")

  object UnknownBackend extends Backend {
    type BackendCall = LocalBackendCall

    def engineFunctions(ioInterface: IoInterface, workflowContext: WorkflowContext): WorkflowEngineFunctions = throw new NotImplementedError

    override val actorSystem = workflowManagerSystem.actorSystem

    override def adjustInputPaths(backendCall: BackendCall) = throw new NotImplementedError

    override def adjustOutputPaths(call: Call, outputs: CallOutputs): CallOutputs = throw new NotImplementedError
    override def stdoutStderr(backendCall: BackendCall): CallLogs = throw new NotImplementedError
    override def initializeForWorkflow(workflow: WorkflowDescriptor) = throw new NotImplementedError
    override def prepareForRestart(restartableWorkflow: WorkflowDescriptor)(implicit ec: ExecutionContext) = throw new NotImplementedError

    override def bindCall(workflowDescriptor: WorkflowDescriptor,
                          key: BackendCallKey,
                          locallyQualifiedInputs: CallInputs,
                          abortRegistrationFunction: Option[AbortRegistrationFunction]): BackendCall =
      throw new NotImplementedError

    override def backendType: BackendType = throw new NotImplementedError
    override def rootPath(workflowOptions: WorkflowOptions): String = throw new NotImplementedError
    override def pollBackoff: ExponentialBackOff = throw new NotImplementedError
    override def executionInfoKeys: List[String] = List.empty
  }

  "SlickDataAccess" should "not deadlock" in {
    // Test based on https://github.com/kwark/slick-deadlock/blob/82525fc/src/main/scala/SlickDeadlock.scala
    val databaseConfig = ConfigFactory.parseString(
      s"""
         |db.url = "jdbc:hsqldb:mem:$${slick.uniqueSchema};shutdown=false;hsqldb.tx=mvcc"
         |db.driver = "org.hsqldb.jdbcDriver"
         |db.numThreads = 2
         |driver = "slick.driver.HsqldbDriver$$"
         |""".stripMargin)

    for {
      dataAccess <- new SlickDataAccess(databaseConfig).autoClosed
    } {
      val futures = 1 to 20 map { _ =>
        val workflowId = WorkflowId.randomId()
        val workflowInfo = WorkflowDescriptor(workflowId, testSources)
        for {
          _ <- dataAccess.createWorkflow(workflowInfo, Nil, Nil, localBackend)
          _ <- dataAccess.getWorkflow(workflowInfo.id) map (_.id should be(workflowId))
        } yield ()
      }
      Future.sequence(futures).futureValue(Timeout(10.seconds))
    }
  }

  it should "have the same liquibase and slick schema" in {
    for {
      liquibaseDataAccess <- dataAccessForSchemaManager("liquibase").autoClosed
      slickDataAccess <- dataAccessForSchemaManager("slick").autoClosed
    } {
      val diffResult = LiquibaseSchemaManager.compare(
        slickDataAccess.dataAccess.driver, slickDataAccess.database,
        liquibaseDataAccess.dataAccess.driver, liquibaseDataAccess.database)

      import DiffResultFilter._
      val diffFilters = StandardTypeFilters :+ UniqueIndexFilter
      val filteredDiffResult = diffResult.filterLiquibaseObjects.filterChangedObjects(diffFilters)

      filteredDiffResult.getChangedObjects should be(empty)
      filteredDiffResult.getMissingObjects should be(empty)
      filteredDiffResult.getUnexpectedObjects should be(empty)
    }
  }

  def dataAccessForSchemaManager(schemaManager: String): SlickDataAccess = {
    val databaseConfig = ConfigFactory.parseString(
      s"""
         |db.url = "jdbc:hsqldb:mem:$${slick.uniqueSchema};shutdown=false;hsqldb.tx=mvcc"
         |db.driver = "org.hsqldb.jdbcDriver"
         |driver = "slick.driver.HsqldbDriver$$"
         |schema.manager = $schemaManager
         |""".stripMargin)
    new SlickDataAccess(databaseConfig)
  }

  // Tests against main database used for command line
  "SlickDataAccess (main.hsqldb)" should behave like databaseWithConfig("main.hsqldb", testRequired = true)

  // Tests using liquibase, but in memory
  "SlickDataAccess (test.hsqldb)" should behave like databaseWithConfig("test.hsqldb", testRequired = true)

  // If able to connect, then also run the tests on mysql, but it's not required
  "SlickDataAccess (test.mysql)" should behave like databaseWithConfig("test.mysql", testRequired = false)

  def databaseWithConfig(path: => String, testRequired: => Boolean): Unit = {

    lazy val testDatabase = new TestSlickDatabase(path)
    lazy val canConnect = testRequired || testDatabase.isValidConnection.futureValue(Timeout(10.seconds))
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

    val testWdlString = WdlString("testStringvalue")
    val testWdlStringHash = testWdlString.computeHash
    val testWdlStringShard = WdlString("testStringValueShard")
    val testWdlStringShardHash = testWdlStringShard.computeHash

    it should "(if hsqldb) have transaction isolation mvcc" taggedAs IntegrationTest in {
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

    it should "create and retrieve the workflow for just reading" taggedAs IntegrationTest in {
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

    it should "store and retrieve an empty String as a WdlValue" taggedAs IntegrationTest in {
      assume(canConnect || testRequired)
      val workflowId = WorkflowId(UUID.randomUUID())
      val workflowInfo = WorkflowDescriptor(workflowId, testSources)
      val call = mock[Call]
      call.fullyQualifiedName returns "wf.a"
      val outputKey = BackendCallKey(call, None, 1)
      val dbKey = ExecutionDatabaseKey(outputKey.scope.fullyQualifiedName, outputKey.index, outputKey.attempt)
      val outputs: CallOutputs = Map("wf.a.empty" -> CallOutput(WdlString(""), None))

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, Nil, Nil, localBackend)
        _ <- dataAccess.setOutputs(workflowId, outputKey, outputs, Seq.empty)
        _ <- dataAccess.getOutputs(workflowId, dbKey) map { results =>
          results shouldNot be(empty)

          val workflowOutput = results.headOption
          workflowOutput shouldNot be(empty)
          val workflowResult = workflowOutput.get
          workflowResult.wdlValue shouldNot be(empty)
          val wdlValue = workflowResult.wdlValue
          wdlValue.get.valueString shouldBe ""
        }
      } yield ()).futureValue
    }

    it should "create and query a workflow" taggedAs IntegrationTest in {
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

    def assertCallCachingFailure(id: WorkflowId, callName: Option[String], messages: String*): Future[Unit] = {
      // The `from` Future is expected to fail, so if the forcomp actually runs the test should fail.
      val parameters = for {
        s <- CallCachingParameters.from(id, None, AllowFalse, dataAccess)
      } yield throw new RuntimeException(s"Unexpected success: $s")

      // `recover` the failed Future looking for an expected `IllegalArgumentException`.  Assert all the expected
      // messages are present in the exception text and the correct number of expected failures are seen.
      // If the `parameters` Future is failed but the exception isn't an `IllegalArgumentException` then this recover
      // won't match and the Future will remain failed and fail the test.
      parameters recover {
        case e: IllegalArgumentException =>
          messages foreach { m => if (!e.getMessage.contains(m)) throw new RuntimeException(s"Missing message: $m.  Exception text: ${e.getMessage}") }
          if (e.getMessage.count(_ == '\n') != messages.size - 1) throw new RuntimeException(s"Unexpected messages seen: ${e.getMessage}")
      }
    }

    it should "support call caching configuration for specified calls in a regular workflow" taggedAs IntegrationTest in {
      assume(canConnect || testRequired)
      val workflowInfo = WorkflowDescriptor(WorkflowId(UUID.randomUUID()), SampleWdl.ThreeStep.asWorkflowSources())

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, Nil, workflowInfo.namespace.workflow.calls, localBackend)
        // Unknown workflow
        _ <- assertCallCachingFailure(WorkflowId(UUID.randomUUID()), callName = Option("three_step.ps"), "Workflow not found")
        _ <- dataAccess.setTerminalStatus(workflowInfo.id, ExecutionDatabaseKey("three_step.ps", None, 1), ExecutionStatus.Done, None, None, None)
        executions <- dataAccess.getExecutions(workflowInfo.id)
        _ = executions should have size 3
        _ = executions foreach { _.allowsResultReuse shouldBe true }
        params <- CallCachingParameters.from(workflowInfo.id, Option("three_step.ps"), AllowFalse, dataAccess)
        _ = params.workflowId shouldBe workflowInfo.id
        _ = params.callKey shouldEqual Option(ExecutionDatabaseKey("three_step.ps", None, 1))
        _ <- dataAccess.updateCallCaching(params)
        executions <- dataAccess.getExecutions(workflowInfo.id)
        (allowing, disallowing) = executions partition { _.allowsResultReuse }
        _ = allowing should have size 2
        _ = disallowing should have size 1
        _ = disallowing.seq.head.callFqn should be("three_step.ps")
      } yield ()).futureValue
    }

    it should "support call caching configuration for specified calls in a scattered workflow" taggedAs IntegrationTest in {
      assume(canConnect || testRequired)
      val workflowInfo = WorkflowDescriptor(WorkflowId(UUID.randomUUID()), SampleWdl.SimpleScatterWdl.asWorkflowSources())

      (for {
        // The `inside_scatter` is a to-be-exploded placeholder, but it will conflict with the collector that the
        // scatter explodes below so filter that out.
        _ <- dataAccess.createWorkflow(workflowInfo, Nil, workflowInfo.namespace.workflow.calls.filterNot(_.unqualifiedName == "inside_scatter"), localBackend)
        scatter = workflowInfo.namespace.workflow.scatters.head
        scatterKey = ScatterKey(scatter, None)
        newEntries = scatterKey.populate(5)
        _ <- dataAccess.insertCalls(workflowInfo.id, newEntries.keys, localBackend)
        executions <- dataAccess.getExecutions(workflowInfo.id)
        _ = executions foreach { _.allowsResultReuse shouldBe true }

        // Calls outside the scatter should work the same as an unscattered workflow.
        outsideParams <- CallCachingParameters.from(workflowInfo.id, Option("scatter0.outside_scatter"), AllowFalse, dataAccess)
        _ <- dataAccess.updateCallCaching(outsideParams)
        executions <- dataAccess.getExecutions(workflowInfo.id)
        (allowing, disallowing) = executions partition { _.allowsResultReuse }
        _ = allowing should have size (executions.size - 1)
        _ = disallowing should have size 1
        _ = disallowing.seq.head.callFqn should be ("scatter0.outside_scatter")

        // Support unindexed scattered call targets to update all shards.
        _ = executions filter { _.callFqn == "scatter0.inside_scatter" } foreach { _.allowsResultReuse shouldBe true }
        unindexedCallParams <- CallCachingParameters.from(workflowInfo.id, Option("scatter0.inside_scatter"), AllowFalse, dataAccess)
        _ <- dataAccess.updateCallCaching(unindexedCallParams)
        executions <- dataAccess.getExecutions(workflowInfo.id)
        _ = executions filter { _.callFqn == "scatter0.inside_scatter" } foreach { _.allowsResultReuse shouldBe false }

        // Support indexed shards as well.
        insideParams <- CallCachingParameters.from(workflowInfo.id, Option("scatter0.inside_scatter.3"), AllowTrue, dataAccess)
        _ <- dataAccess.updateCallCaching(insideParams)
        executions <- dataAccess.getExecutions(workflowInfo.id)
        inside = executions filter { e => e.callFqn == "scatter0.inside_scatter" }
        (allowing, disallowing) = inside partition { _.allowsResultReuse }
        _ = allowing should have size 1
        _ = disallowing should have size (inside.size - 1)
      } yield ()).futureValue
    }

    it should "support call caching configuration for all calls in a regular workflow" taggedAs IntegrationTest in {
      assume(canConnect || testRequired)

      val workflowInfo = WorkflowDescriptor(WorkflowId(UUID.randomUUID()), SampleWdl.ThreeStep.asWorkflowSources())
      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, Nil, workflowInfo.namespace.workflow.calls, localBackend)
        // Unknown workflow
        _ <- assertCallCachingFailure(WorkflowId(UUID.randomUUID()), callName = None, "Workflow not found")
        params <- CallCachingParameters.from(workflowInfo.id, None, AllowFalse, dataAccess)
        executions <- dataAccess.getExecutions(workflowInfo.id)
        _ = executions should have size 3
        _ = executions foreach { _.allowsResultReuse shouldBe true }
        _ <- dataAccess.updateCallCaching(params)
        executions <- dataAccess.getExecutions(workflowInfo.id)
        _ = executions foreach { _.allowsResultReuse shouldBe false }
      } yield ()).futureValue
    }

    it should "support call caching configuration for all calls in a scattered workflow" taggedAs IntegrationTest in {
      assume(canConnect || testRequired)

      val workflowInfo = WorkflowDescriptor(WorkflowId(UUID.randomUUID()), SampleWdl.SimpleScatterWdl.asWorkflowSources())
      (for {
        // The `inside_scatter` is a to-be-exploded placeholder, but it will conflict with the collector that the
        // scatter explodes below so filter that out.
        _ <- dataAccess.createWorkflow(workflowInfo, Nil, workflowInfo.namespace.workflow.calls.filterNot(_.unqualifiedName == "inside_scatter"), localBackend)
        scatter = workflowInfo.namespace.workflow.scatters.head
        scatterKey = ScatterKey(scatter, None)
        newEntries = scatterKey.populate(5)
        _ <- dataAccess.insertCalls(workflowInfo.id, newEntries.keys, localBackend)

        scatterParams <- CallCachingParameters.from(workflowInfo.id, None, AllowFalse, dataAccess)
        executions <- dataAccess.getExecutions(workflowInfo.id)
        _ = executions foreach { _.allowsResultReuse shouldBe true }
        _ <- dataAccess.updateCallCaching(scatterParams)
        executions <- dataAccess.getExecutions(workflowInfo.id)
        _ = executions foreach { _.allowsResultReuse shouldBe false }
      } yield ()).futureValue
    }


    it should "query a single execution status" taggedAs IntegrationTest in {
      assume(canConnect || testRequired)
      val workflowId = WorkflowId(UUID.randomUUID())
      val workflowInfo = WorkflowDescriptor(workflowId, testSources)
      val task = Task.empty
      val callFqn = "fully.qualified.name"
      val call = new Call(None, callFqn, task, Set.empty[FullyQualifiedName], Map.empty, None)

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, Nil, Seq(call), localBackend)
        _ <- dataAccess.updateWorkflowState(workflowId, WorkflowRunning)
        _ <- dataAccess.getExecutionStatus(workflowId, ExecutionDatabaseKey(callFqn, None, 1)) map { status =>
          status.get.executionStatus shouldBe ExecutionStatus.NotStarted
          status.get.returnCode shouldBe None
        }
      } yield ()).futureValue
    }

    it should "create and retrieve 3step.wdl with a 10,000 char pattern" taggedAs IntegrationTest in {
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

    it should "fail when saving a workflow twice" taggedAs IntegrationTest in {
      assume(canConnect || testRequired)
      val workflowId = WorkflowId(UUID.randomUUID())
      val workflowInfo = WorkflowDescriptor(workflowId, testSources)

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, Nil, Nil, localBackend)
        _ <- dataAccess.createWorkflow(workflowInfo, Nil, Nil, localBackend)
      } yield ()).failed.futureValue should be(a[SQLException])
    }

    it should "fail when updating a non-existent workflow state" taggedAs IntegrationTest in {
      assume(canConnect || testRequired)
      val workflowId = WorkflowId(UUID.randomUUID())

      (for {
        _ <- dataAccess.updateWorkflowState(workflowId, WorkflowRunning)
      } yield ()).failed.futureValue should be(an[IllegalArgumentException])
    }

    it should "update and get a workflow state" taggedAs IntegrationTest in {
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

    it should "get workflow state" taggedAs IntegrationTest in {
      assume(canConnect || testRequired)
      val workflowId = WorkflowId(UUID.randomUUID())
      val workflowInfo = WorkflowDescriptor(workflowId, testSources)
      val key = new SymbolStoreKey("myScope", "myName", None, input = true)
      val entry = new SymbolStoreEntry(key, WdlStringType, Option(testWdlString), Option(testWdlStringHash))

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

      it should s"get $spec" taggedAs IntegrationTest in {
        assume(canConnect || testRequired)
        val workflowId = WorkflowId(UUID.randomUUID())
        val workflowInfo = WorkflowDescriptor(workflowId, testSources)

        val task = Task.empty
        val call =
          if (setCallParent) {
            val workflow = new Workflow("workflow.name", Nil, Nil)
            val call = new Call(callAlias, "call.name", task, Set.empty[FullyQualifiedName], Map.empty, Option(workflow))
            workflow.children = Seq(call)
            call
          } else {
            new Call(callAlias, "call.name", task, Set.empty[FullyQualifiedName], Map.empty, None)
          }

        val callKey = ExecutionDatabaseKey(call.fullyQualifiedName, None, 1)
        val shardKey = ExecutionDatabaseKey(call.fullyQualifiedName, Option(0), 1)

        def optionallyUpdateExecutionStatus() =
          if (updateStatus)
            dataAccess.updateStatus(workflowId, Seq(callKey), ExecutionStatus.Running)
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
          _ <- dataAccess.insertCalls(workflowId, Seq(BackendCallKey(call, Option(0), 1)), localBackend)
          _ <- dataAccess.setTerminalStatus(workflowId, shardKey, ExecutionStatus.Done, Option(0), None, None)
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

    it should "insert a call in execution table" taggedAs IntegrationTest in {
      assume(canConnect || testRequired)
      val callFqn = "call.fully.qualified.scope"
      val workflowId = WorkflowId(UUID.randomUUID())
      val workflowInfo = WorkflowDescriptor(workflowId, testSources)
      val task = Task.empty
      val call = new Call(None, callFqn, task, Set.empty[FullyQualifiedName], Map.empty, None)
      val callKey = BackendCallKey(call, None, 1)
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

    it should "fail to get an non-existent execution status" taggedAs IntegrationTest in {
      assume(canConnect || testRequired)
      dataAccess.getExecutionStatuses(WorkflowId(UUID.randomUUID())).failed.futureValue should be(a[NoSuchElementException])
    }

    it should "get a symbol input" taggedAs IntegrationTest in {
      assume(canConnect || testRequired)
      val callFqn = "call.fully.qualified.scope"
      val symbolFqn = "symbol.fully.qualified.scope"
      val workflowId = WorkflowId(UUID.randomUUID())
      val workflowInfo = WorkflowDescriptor(workflowId, testSources)
      val key = new SymbolStoreKey(callFqn, symbolFqn, None, input = true)
      val entry = new SymbolStoreEntry(key, WdlStringType, Option(testWdlString), Option(testWdlStringHash))
      val task = Task.empty
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
          resultSymbol.wdlValue.get should be(testWdlString)
          resultSymbol.symbolHash shouldNot be(empty)
          resultSymbol.symbolHash should be(Option(testWdlStringHash))
        }
      } yield ()).futureValue
    }

    it should "get a symbol input that has a very long WDL value field" taggedAs IntegrationTest in {
      assume(canConnect || testRequired)
      val wdlArray = new WdlArray(WdlArrayType(WdlStringType), Seq(WdlString("test"), WdlString("*" * 10000)))
      val callFqn = "call.fully.qualified.scope"
      val symbolFqn = "symbol.fully.qualified.scope"
      val workflowId = WorkflowId(UUID.randomUUID())
      val workflowDescriptor = WorkflowDescriptor(workflowId, testSources)
      val key = new SymbolStoreKey(callFqn, symbolFqn, None, input = true)
      val entry = new SymbolStoreEntry(key, WdlArrayType(WdlStringType), Option(wdlArray), wdlArray.getHash(workflowDescriptor))
      val task = Task.empty
      val call = new Call(None, callFqn, task, Set.empty[FullyQualifiedName], Map.empty, None)

      (for {
        _ <- dataAccess.createWorkflow(workflowDescriptor, Seq(entry), Nil, localBackend)
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
          resultSymbol.symbolHash shouldNot be(empty)
          resultSymbol.symbolHash should be(Option(wdlArray.computeHash))
        }
      } yield ()).futureValue
    }

    it should "fail to get inputs for a null call" taggedAs IntegrationTest in {
      assume(canConnect || testRequired)
      val workflowId: WorkflowId = WorkflowId(UUID.randomUUID())
      val workflowInfo = WorkflowDescriptor(workflowId, testSources)

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, Nil, Nil, localBackend)
        _ <- dataAccess.updateWorkflowState(workflowId, WorkflowRunning)
        _ <- dataAccess.getInputs(workflowId, null)
      } yield ()).failed.futureValue should be(an[IllegalArgumentException])
    }

    it should "set and get an output" taggedAs IntegrationTest in {
      assume(canConnect || testRequired)
      val callFqn = "call.fully.qualified.scope"
      val symbolLqn = "symbol"
      val workflowId = WorkflowId(UUID.randomUUID())
      val workflowInfo = WorkflowDescriptor(workflowId, testSources)
      val task = Task.empty
      val call = new Call(None, callFqn, task, Set.empty[FullyQualifiedName], Map.empty, None)

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, Nil, Nil, localBackend)
        _ <- dataAccess.updateWorkflowState(workflowId, WorkflowRunning)
        _ <- dataAccess.setOutputs(workflowId, BackendCallKey(call, None, 1), Map(symbolLqn -> CallOutput(testWdlString, Option(testWdlStringHash))), Seq.empty)
        _ <- dataAccess.setOutputs(workflowId, BackendCallKey(call, Option(0), 1), Map(symbolLqn -> CallOutput(testWdlStringShard, Option(testWdlStringShardHash))), Seq.empty)
        _ <- dataAccess.getOutputs(workflowId, ExecutionDatabaseKey(call.fullyQualifiedName, None, 1)) map { results =>
          results.size should be(1) //getOutputs on a workflowId does NOT return shards outputs

          val output = results.head
          output.key.scope should be("call.fully.qualified.scope")
          output.key.name should be("symbol")
          output.key.input should be(right = false) // IntelliJ highlighting
          output.wdlType should be(WdlStringType)
          output.wdlValue shouldNot be(empty)
          output.key.index should be(None)
          output.wdlValue.get should be(testWdlString)
        }
      } yield ()).futureValue
    }

    it should "set and get shard statuses" taggedAs IntegrationTest in {
      assume(canConnect || testRequired)
      val callFqn1 = "call.fully.qualified.scope$s1"
      val callFqn2 = "call.fully.qualified.scope$s2"
      val workflowId = WorkflowId(UUID.randomUUID())
      val workflowInfo = WorkflowDescriptor(workflowId, testSources)
      val task = Task.empty
      val call1 = new Call(None, callFqn1, task, Set.empty[FullyQualifiedName], Map.empty, None)
      val shardIndex1 = Option(0)
      val pid1 = Option("123")
      val call2 = new Call(None, callFqn2, task, Set.empty[FullyQualifiedName], Map.empty, None)
      val shardIndex2 = Option(1)
      val pid2 = Option("987")

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, Nil, Nil, localBackend)
        _ <- dataAccess.insertCalls(workflowId, Seq(BackendCallKey(call1, shardIndex1, 1)), localBackend)
        _ <- dataAccess.insertCalls(workflowId, Seq(BackendCallKey(call2, shardIndex2, 1)), localBackend)
        _ <- dataAccess.updateExecutionInfo(workflowId, BackendCallKey(call1, shardIndex1, 1), InfoKeys.Pid, pid1)
        _ <- dataAccess.updateExecutionInfo(workflowId, BackendCallKey(call2, shardIndex2, 1), InfoKeys.Pid, pid2)
        _ <- dataAccess.getExecutionInfos(workflowId, call1, 1) map {
          case infos => assertResult(pid1.get) { infos.head.value.get }
        }
        _ <- dataAccess.getExecutionInfos(workflowId, call2, 1) map {
          case infos => assertResult(pid2.get) { infos.head.value.get }
        }
      } yield ()).futureValue
    }

    it should "set and get an output by call" taggedAs IntegrationTest in {
      assume(canConnect || testRequired)
      val callFqn = "call.fully.qualified.scope"
      val symbolLqn = "symbol"
      val workflowId = WorkflowId(UUID.randomUUID())
      val workflowInfo = WorkflowDescriptor(workflowId, testSources)
      val task = Task.empty
      val call = new Call(None, callFqn, task, Set.empty[FullyQualifiedName], Map.empty, None)

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, Nil, Nil, localBackend)
        _ <- dataAccess.updateWorkflowState(workflowId, WorkflowRunning)
        _ <- dataAccess.setOutputs(workflowId, BackendCallKey(call, None, 1), Map(symbolLqn -> CallOutput(testWdlString, Option(testWdlStringHash))), Seq.empty)
        _ <- dataAccess.setOutputs(workflowId, BackendCallKey(call, Option(0), 1), Map(symbolLqn -> CallOutput(testWdlStringShard, Option(testWdlStringShardHash))), Seq.empty)
        callOutput <- dataAccess.getOutputs(workflowId, ExecutionDatabaseKey(call.fullyQualifiedName, None, 1)) map { results =>
          results.head.key.index should be(None)
          results.head.wdlValue.get should be(testWdlString)
          results
        }
        shardOutput <- dataAccess.getOutputs(workflowId, ExecutionDatabaseKey(call.fullyQualifiedName, Option(0), 1)) map { results =>
          results.head.key.index should be(Some(0))
          results.head.wdlValue.get should be(testWdlStringShard)
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

    it should "fail to get outputs for a null call" taggedAs IntegrationTest in {
      assume(canConnect || testRequired)
      val workflowId = WorkflowId(UUID.randomUUID())
      val workflowInfo = WorkflowDescriptor(workflowId, testSources)

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, Nil, Nil, localBackend)
        _ <- dataAccess.updateWorkflowState(workflowId, WorkflowRunning)
        _ <- dataAccess.getOutputs(workflowId, null)
      } yield ()).failed.futureValue should be(an[IllegalArgumentException])
    }

    it should "fail to create workflow for an unknown backend" taggedAs IntegrationTest in {
      assume(canConnect || testRequired)
      val callFqn = "call.fully.qualified.scope"
      val workflowId = WorkflowId(UUID.randomUUID())
      val workflowInfo = WorkflowDescriptor(workflowId, testSources)
      val task = Task.empty
      val call = new Call(None, callFqn, task, Set.empty[FullyQualifiedName], Map.empty, None)

      dataAccess.createWorkflow(workflowInfo, Nil, Seq(call),
        UnknownBackend).failed.futureValue should be(an[Exception])
    }

    it should "fail to create workflow for a null backend" taggedAs IntegrationTest in {
      assume(canConnect || testRequired)
      val callFqn = "call.fully.qualified.scope"
      val workflowId = WorkflowId(UUID.randomUUID())
      val workflowInfo = WorkflowDescriptor(workflowId, testSources)
      val task = Task.empty
      val call = new Call(None, callFqn, task, Set.empty[FullyQualifiedName], Map.empty, None)

      dataAccess.createWorkflow(workflowInfo, Nil, Seq(call),
        null).failed.futureValue should be(a[NullPointerException])
    }

    it should "set and get the same symbol with IO as input then output" taggedAs IntegrationTest in {
      assume(canConnect || testRequired)
      val callFqn = "call.fully.qualified.scope"
      val symbolLqn = "symbol"
      val symbolFqn = callFqn + "." + symbolLqn
      val workflowId = WorkflowId(UUID.randomUUID())
      val workflowInfo = WorkflowDescriptor(workflowId, testSources)
      val key = new SymbolStoreKey(callFqn, symbolFqn, None, input = true)
      val entry = new SymbolStoreEntry(key, WdlStringType, Option(testWdlString), Option(testWdlStringHash))
      val task = Task.empty
      val call = new Call(None, callFqn, task, Set.empty[FullyQualifiedName], Map.empty, None)

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, Seq(entry), Nil, localBackend)
        _ <- dataAccess.updateWorkflowState(workflowId, WorkflowRunning)
        _ <- dataAccess.setOutputs(workflowId, BackendCallKey(call, None, 1), Map(symbolLqn -> CallOutput(testWdlString, Option(testWdlStringHash))), Seq.empty)
        _ <- dataAccess.getOutputs(workflowId, ExecutionDatabaseKey(call.fullyQualifiedName, None, 1)) map { results =>
          results.size should be(1)
          val resultSymbol = results.head
          val resultSymbolStoreKey = resultSymbol.key
          resultSymbolStoreKey.scope should be("call.fully.qualified.scope")
          resultSymbolStoreKey.name should be("symbol")
          resultSymbolStoreKey.index should be(None)
          resultSymbolStoreKey.input should be(right = false) // IntelliJ highlighting
          resultSymbol.wdlType should be(WdlStringType)
          resultSymbol.wdlValue shouldNot be(empty)
          resultSymbol.wdlValue.get should be(testWdlString)
          resultSymbol.symbolHash shouldNot be(empty)
          resultSymbol.symbolHash should be(Option(testWdlStringHash))
        }
      } yield ()).futureValue
    }

    it should "fail when setting an existing symbol output" taggedAs IntegrationTest in {
      assume(canConnect || testRequired)
      val callFqn = "call.fully.qualified.scope"
      val symbolLqn = "symbol"
      val symbolFqn = callFqn + "." + symbolLqn
      val workflowId = WorkflowId(UUID.randomUUID())
      val workflowInfo = WorkflowDescriptor(workflowId, testSources)
      val task = Task.empty
      val call = new Call(None, callFqn, task, Set.empty[FullyQualifiedName], Map.empty, None)

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, Seq(), Nil, localBackend)
        _ <- dataAccess.updateWorkflowState(workflowId, WorkflowRunning)
        _ <- dataAccess.setOutputs(workflowId, BackendCallKey(call, None, 1), Map(symbolFqn -> CallOutput(testWdlString, Option(testWdlStringHash))), Seq.empty)
        // Second attempt should fail
        _ <- dataAccess.setOutputs(workflowId, BackendCallKey(call, None, 1), Map(symbolFqn -> CallOutput(testWdlString, Option(testWdlStringHash))), Seq.empty)
      } yield ()).failed.futureValue should be(a[SQLException])
    }

    it should "set and get a backend info" taggedAs IntegrationTest in {
      assume(canConnect || testRequired)
      val workflowId = WorkflowId(UUID.randomUUID())
      val workflowInfo = WorkflowDescriptor(workflowId, testSources)
      val task = Task.empty
      val call = new Call(None, "fully.qualified.name", task, Set.empty[FullyQualifiedName], Map.empty, None)
      val pid = Option("123")

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, Nil, Seq(call), localBackend)
        _ <- dataAccess.updateExecutionInfo(workflowId, BackendCallKey(call, None, 1), InfoKeys.Pid, pid)
        _ <- dataAccess.getExecutionInfos(workflowId, call, 1) map { info =>
          info should have size 1
          info.head.key should equal(InfoKeys.Pid)
          info.head.value should equal(pid)
        }
      } yield ()).futureValue
    }

    // Queries use `.head` a lot. There was a bug that pulled the backend info by fqn, but for any workflow.
    it should "set and get a backend info for same call on two workflows" taggedAs IntegrationTest in {
      assume(canConnect || testRequired)
      val workflowId1 = WorkflowId(UUID.randomUUID())
      val workflowId2 = WorkflowId(UUID.randomUUID())
      val workflowInfo1 = WorkflowDescriptor(workflowId1, testSources)
      val workflowInfo2 = WorkflowDescriptor(workflowId2, testSources)
      val task = Task.empty
      val call = new Call(None, "fully.qualified.name", task, Set.empty[FullyQualifiedName], Map.empty, None)
      val pid1 = Option("123")
      val pid2 = Option("321")

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo1, Nil, Seq(call), localBackend)
        _ <- dataAccess.createWorkflow(workflowInfo2, Nil, Seq(call), localBackend)
        _ <- dataAccess.updateExecutionInfo(workflowId1, BackendCallKey(call, None, 1), InfoKeys.Pid, pid1)
        _ <- dataAccess.updateExecutionInfo(workflowId2, BackendCallKey(call, None, 1), InfoKeys.Pid, pid2)
        _ <- dataAccess.getExecutionInfos(workflowId1, call, 1) map { info =>
          info should have size 1
          info.head.key should equal(InfoKeys.Pid)
          info.head.value should equal(pid1)
        }
        _ <- dataAccess.getExecutionInfos(workflowId2, call, 1) map { info =>
          info should have size 1
          info.head.key should equal(InfoKeys.Pid)
          info.head.value should equal(pid2)
        }
      } yield ()).futureValue
    }

    it should "update call start and end dates appropriately" taggedAs IntegrationTest in {
      assume(canConnect || testRequired)
      val callFqn = "call.fully.qualified.scope"
      val workflowId = WorkflowId(UUID.randomUUID())
      val workflowInfo = WorkflowDescriptor(workflowId, testSources)
      val task = Task.empty
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
        callKey = ExecutionDatabaseKey(callFqn, None, 1)

        _ <- dataAccess.updateStatus(workflowId, List(callKey), ExecutionStatus.NotStarted)
        _ <- dataAccess.getExecutions(workflowId) map assertDates(startDefined = false, endDefined = false)

        _ <- dataAccess.setStartingStatus(workflowId, List(callKey))
        _ <- dataAccess.getExecutions(workflowId) map assertDates(startDefined = true, endDefined = false)

        _ <- dataAccess.setTerminalStatus(workflowId, callKey, ExecutionStatus.Done, None, None, None)
        _ <- dataAccess.getExecutions(workflowId) map assertDates(startDefined = true, endDefined = true)
      } yield()).futureValue
    }

    it should "set and get execution events" taggedAs IntegrationTest in {
      assume(canConnect || testRequired)

      // We need an execution to create an event. We need a workflow to make an execution. Le Sigh...
      val workflowId = WorkflowId(UUID.randomUUID())
      val workflowInfo = WorkflowDescriptor(workflowId, testSources)
      val task = Task.empty
      val call = new Call(None, "fully.qualified.name", task, Set.empty[FullyQualifiedName], Map.empty, None)
      val shardedCall = new Call(None, "fully.qualified.name", task, Set.empty[FullyQualifiedName], Map.empty, None)
      val shardIndex = Some(0)

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, Nil, Seq(call), localBackend)
        _ <- dataAccess.insertCalls(workflowId, Seq(BackendCallKey(call, None, 2)), localBackend)
        _ <- dataAccess.insertCalls(workflowId, Seq(BackendCallKey(shardedCall, shardIndex, 1)), localBackend)
        _ <- dataAccess.insertCalls(workflowId, Seq(BackendCallKey(shardedCall, shardIndex, 2)), localBackend)

        now = DateTime.now
        mainEventSeq = Seq(
          ExecutionEventEntry("hello", now.minusHours(7), now.minusHours(5)),
          ExecutionEventEntry("o-genki desuka", now.minusHours(4), now.minusHours(2)),
          ExecutionEventEntry("oopsPreempted", now.minusHours(5), now))
        mainEventSeqAttempt2 = Seq(
          ExecutionEventEntry("hello", now.minusHours(7), now.minusHours(5)),
          ExecutionEventEntry("cheerio", now.minusHours(5), now))
        shardEventSeq = Seq(
          ExecutionEventEntry("konnichiwa", now.minusHours(4), now.minusHours(2)),
          ExecutionEventEntry("preempted again!", now.minusHours(1), now))
        shardEventSeqAttempt2 = Seq(
          ExecutionEventEntry("konnichiwa", now.minusHours(7), now.minusHours(5)),
          ExecutionEventEntry("o-genki desuka", now.minusHours(5), now.minusHours(2)),
          ExecutionEventEntry("sayounara", now.minusHours(2), now))

        _ <- dataAccess.setExecutionEvents(workflowId, call.fullyQualifiedName, None, 1, mainEventSeq)
        _ <- dataAccess.setExecutionEvents(workflowId, call.fullyQualifiedName, None, 2, mainEventSeqAttempt2)
        _ <- dataAccess.setExecutionEvents(workflowId, call.fullyQualifiedName, shardIndex, 1, shardEventSeq)
        _ <- dataAccess.setExecutionEvents(workflowId, call.fullyQualifiedName, shardIndex, 2, shardEventSeqAttempt2)

        _ <- dataAccess.getAllExecutionEvents(workflowId) map { retrievedEvents =>
          val mainExecutionDatabaseKey = ExecutionDatabaseKey(call.fullyQualifiedName, None, 1)
          retrievedEvents valueAt mainExecutionDatabaseKey should have size 3
          mainEventSeq foreach { event =>
            retrievedEvents valueAt mainExecutionDatabaseKey exists { executionEventsCloseEnough(_, event) } should be (true)
          }

          val mainExecutionDatabaseKeyAttempt2 = ExecutionDatabaseKey(call.fullyQualifiedName, None, 2)
          retrievedEvents valueAt mainExecutionDatabaseKeyAttempt2 should have size 2
          mainEventSeqAttempt2 foreach { event =>
            retrievedEvents valueAt mainExecutionDatabaseKeyAttempt2 exists { executionEventsCloseEnough(_, event) } should be (true)
          }

          val shardExecutionDatabaseKey = ExecutionDatabaseKey(call.fullyQualifiedName, shardIndex, 1)
          retrievedEvents valueAt shardExecutionDatabaseKey should have size 2
          shardEventSeq foreach { event =>
            retrievedEvents valueAt shardExecutionDatabaseKey exists { executionEventsCloseEnough(_, event) } should be (true)
          }

          val shardExecutionDatabaseKeyAttempt2 = ExecutionDatabaseKey(call.fullyQualifiedName, shardIndex, 2)
          retrievedEvents valueAt shardExecutionDatabaseKeyAttempt2 should have size 3
          shardEventSeqAttempt2 foreach { event =>
            retrievedEvents valueAt shardExecutionDatabaseKeyAttempt2 exists { executionEventsCloseEnough(_, event) } should be (true)
          }
        }
      } yield ()).futureValue
    }

    it should "reject a set of execution events without a valid execution to link to" taggedAs IntegrationTest in {
      assume(canConnect || testRequired)

      // We need an execution to create an event. We need a workflow to make an execution. Le Sigh...
      val workflowId = WorkflowId(UUID.randomUUID())
      val workflowInfo = WorkflowDescriptor(workflowId, testSources)
      val task = Task.empty
      val call = new Call(None, "fully.qualified.name", task, Set.empty[FullyQualifiedName], Map.empty, None)

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, Nil, Seq(call), localBackend)
        _ <- dataAccess.setExecutionEvents(WorkflowId(UUID.randomUUID()), call.fullyQualifiedName, None, 1, Seq(
          ExecutionEventEntry("hello", DateTime.now.minusHours(7), DateTime.now.minusHours(5)),
          ExecutionEventEntry("cheerio", DateTime.now.minusHours(5), DateTime.now)))
      } yield ()).failed.futureValue should be(a[NoSuchElementException])
    }

    it should "close the database" taggedAs IntegrationTest in {
      assume(canConnect || testRequired)
      dataAccess.close()
    }
  }

  private def executionEventsCloseEnough(a: ExecutionEventEntry, b: ExecutionEventEntry): Boolean = {
    a.description == b.description &&
      a.startTime.getMillis - a.startTime.getMillisOfSecond == b.startTime.getMillis - b.startTime.getMillisOfSecond &&
      a.endTime.getMillis - a.endTime.getMillisOfSecond == b.endTime.getMillis - b.endTime.getMillisOfSecond
  }
}
