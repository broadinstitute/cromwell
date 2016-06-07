package cromwell.engine.db.slick

import java.nio.file.{FileSystem, FileSystems}
import java.sql.SQLException
import java.time.OffsetDateTime

import better.files._
import com.google.api.client.util.ExponentialBackOff
import com.typesafe.config.ConfigFactory
import cromwell.CromwellSpec.DbmsTest
import cromwell.CromwellTestkitSpec.TestWorkflowManagerSystem
import cromwell.backend.wdl.{OldCallEngineFunctions, OldWorkflowEngineFunctions}
import cromwell.backend.JobKey
import cromwell.core._
import cromwell.database.SqlConverters._
import cromwell.database.obj.{Execution, WorkflowMetadataKeys}
import cromwell.database.slick.SlickDatabase
import cromwell.engine._
import cromwell.engine.backend._
import cromwell.engine.backend.io._
import cromwell.engine.backend.local.OldStyleLocalBackend
import cromwell.engine.backend.local.OldStyleLocalBackend.InfoKeys
import cromwell.engine.db.{DataAccess, ExecutionDatabaseKey}
import cromwell.engine.workflow.{BackendCallKey, ScatterKey}
import cromwell.services.{MetadataEvent, MetadataKey, MetadataValue}
import cromwell.util.SampleWdl
import cromwell.webservice.{CallCachingParameters, WorkflowQueryKey, WorkflowQueryParameters}
import cromwell.{CromwellTestkitSpec, webservice}
import org.scalactic.StringNormalizations._
import org.scalatest.PartialFunctionValues._
import org.scalatest.concurrent.PatienceConfiguration.Timeout
import org.scalatest.concurrent.ScalaFutures
import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.time.SpanSugar._
import org.scalatest.time.{Millis, Seconds, Span}
import org.scalatest.{BeforeAndAfterAll, FlatSpec, Matchers}
import org.specs2.mock.Mockito
import wdl4s._
import wdl4s.types.{WdlArrayType, WdlStringType}
import wdl4s.values.{WdlFile, _}

import scala.concurrent.{ExecutionContext, Future}
import scala.util.Try

object SlickDataAccessSpec {
  val AllowFalse = Seq(webservice.QueryParameter("allow", "false"))
  val AllowTrue = Seq(webservice.QueryParameter("allow", "true"))

  val Workflow1Name = "test1"
  val Workflow2Name = "test2"

  implicit class EnhancedEngineWorkflowDescriptor(val workflowDescriptor: EngineWorkflowDescriptor) extends AnyVal {
    def fileHasher: FileHasher = { wdlFile: WdlFile =>
        SymbolHash(wdlFile.value.toPath(FileSystems.getDefault).hash)
    }
  }
}

class SlickDataAccessSpec extends FlatSpec with Matchers with ScalaFutures with BeforeAndAfterAll with Mockito
  with WorkflowDescriptorBuilder {
  import SlickDataAccessSpec._

  behavior of "SlickDataAccess"

  import TableDrivenPropertyChecks._

  val workflowManagerSystem = new TestWorkflowManagerSystem
  override implicit val actorSystem = workflowManagerSystem.actorSystem

  override protected def afterAll() = {
    workflowManagerSystem.shutdownTestActorSystem()
    super.afterAll()
  }

  implicit val ec = ExecutionContext.global

  implicit val defaultPatience = PatienceConfig(timeout = Span(5, Seconds), interval = Span(100, Millis))

  lazy val localBackend = OldStyleLocalBackend(CromwellTestkitSpec.DefaultLocalBackendConfigEntry, workflowManagerSystem.actorSystem)

  val test1Sources = WorkflowSourceFiles("workflow test1 {}", "{}", "{}")
  val test2Sources = WorkflowSourceFiles("workflow test2 {}", "{}", "{}")

  object UnknownOldStyleBackend$ extends OldStyleBackend {
    def engineFunctions(fileSystems: List[FileSystem], workflowContext: OldWorkflowContext): OldWorkflowEngineFunctions = throw new NotImplementedError
    override val actorSystem = workflowManagerSystem.actorSystem
    override val backendConfigEntry = CromwellTestkitSpec.DefaultLocalBackendConfigEntry
    override def adjustInputPaths(backendCallJobDescriptor: OldStyleBackendCallJobDescriptor) = throw new NotImplementedError()
    override def stdoutStderr(jobDescriptor: OldStyleBackendCallJobDescriptor): CallLogs = throw new NotImplementedError
    override def initializeForWorkflow(workflow: OldStyleWorkflowDescriptor) = throw new NotImplementedError
    override def backendType: BackendType = throw new NotImplementedError
    override def rootPath(workflowOptions: WorkflowOptions): String = throw new NotImplementedError
    override def pollBackoff: ExponentialBackOff = throw new NotImplementedError
    override def executionInfoKeys: List[String] = List.empty
    override def instantiateCommand(descriptor: OldStyleBackendCallJobDescriptor): Try[String] = throw new NotImplementedError
    override def poll(jobDescriptor: OldStyleBackendCallJobDescriptor, previous: OldStyleExecutionHandle)(implicit ec: ExecutionContext): Future[OldStyleExecutionHandle] = throw new NotImplementedError()
    override def callEngineFunctions(descriptor: OldStyleBackendCallJobDescriptor): OldCallEngineFunctions = throw new NotImplementedError()
    override def useCachedCall(cachedCall: OldStyleBackendCallJobDescriptor, backendCall: OldStyleBackendCallJobDescriptor)(implicit ec: ExecutionContext): Future[OldStyleExecutionHandle] = throw new NotImplementedError()
    override def execute(jobDescriptor: OldStyleBackendCallJobDescriptor)(implicit ec: ExecutionContext): Future[OldStyleExecutionHandle] = throw new NotImplementedError()
    override def resume(descriptor: OldStyleBackendCallJobDescriptor, executionInfos: Map[String, Option[String]])(implicit ec: ExecutionContext): Future[OldStyleExecutionHandle] = throw new NotImplementedError()
    override def fileSystems(options: WorkflowOptions): List[FileSystem] = List(FileSystems.getDefault)
    override def isResumable(key: JobKey, executionInfo: Map[String, Option[String]]) = false
    override def isRestartable(key: JobKey, executionInfo: Map[String, Option[String]]) = false
  }

  it should "not deadlock" in {
    // Test based on https://github.com/kwark/slick-deadlock/blob/82525fc/src/main/scala/SlickDeadlock.scala
    val databaseConfig = ConfigFactory.parseString(
      s"""
         |db.url = "jdbc:hsqldb:mem:$${slick.uniqueSchema};shutdown=false;hsqldb.tx=mvcc"
         |db.driver = "org.hsqldb.jdbcDriver"
         |db.numThreads = 2
         |driver = "slick.driver.HsqldbDriver$$"
         |""".stripMargin)

    for {
      dataAccess <- (new SlickDatabase(databaseConfig) with DataAccess).autoClosed
    } {
      val futures = 1 to 20 map { _ =>
        val workflowId = WorkflowId.randomId()
        val workflowInfo = createMaterializedEngineWorkflowDescriptor(id = workflowId, workflowSources = test1Sources)
        for {
          _ <- dataAccess.createWorkflow(workflowInfo, test1Sources, Nil, Nil, localBackend)
          _ <- dataAccess.getWorkflowExecutionAndAux(workflowInfo.id) map { result =>
            result.execution.workflowExecutionUuid should be(workflowId.toString)
          }
        } yield ()
      }
      Future.sequence(futures).futureValue(Timeout(10.seconds))
    }
  }

  "SlickDataAccess (main.hsqldb)" should behave like testWith("main.hsqldb")

  "SlickDataAccess (test.mysql)" should behave like testWith("test.mysql")

  /*
  If the above tests fail, one may be accidentally talking to the singleton.
  Uncomment the following test, and see if the suite works. If so, have any code talking to services etc. talk to
  the DataAccess instead.
  */
  //"SlickDataAccess (global.singleton)" should behave like testWith("global.singleton")

  def testWith(configPath: String): Unit = {

    lazy val dataAccess: SlickDatabase with DataAccess = configPath match {
      case "global.singleton" => DataAccess.globalDataAccess.asInstanceOf[SlickDatabase with DataAccess]
      case "main.hsqldb" => new SlickDatabase() with DataAccess // Test the no-args constructor
      case _ => new SlickDatabase(SlickDatabase.getDatabaseConfig(configPath)) with DataAccess
    }

    /**
      * Assert that the specified workflow has the appropriate number of calls in appropriately terminal states per
      * the `ExpectedTerminal` function.  This function maps from an FQN + index pair to a Boolean expectation.
      */
    type ExpectedTerminal = (FullyQualifiedName, Int) => Boolean
    def assertTerminalExecution(id: WorkflowId, expectedCount: Int, expectedTerminal: ExpectedTerminal): Future[Unit] = {
      for {
        _ <- dataAccess.getExecutions(id.toString) map { executions =>
          executions should have size expectedCount
          executions foreach { execution =>
            execution.endDt.isDefined shouldBe expectedTerminal(execution.callFqn, execution.index)
          }
        }
      } yield ()
    }

    val id = WorkflowId.randomId()
    import SlickDataAccessSpec.EnhancedEngineWorkflowDescriptor
    implicit lazy val hasher = createMaterializedEngineWorkflowDescriptor(id, test1Sources).fileHasher
    lazy val testWdlString = WdlString("testStringvalue")
    lazy val testWdlStringHash = testWdlString.computeHash
    lazy val testWdlStringShard = WdlString("testStringValueShard")
    lazy val testWdlStringShardHash = testWdlStringShard.computeHash

    it should "(if hsqldb) have transaction isolation mvcc" taggedAs DbmsTest in {
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

    it should "create and retrieve the workflow for just reading" taggedAs DbmsTest in {
      val workflowId = WorkflowId.randomId()
      val workflowInfo = createMaterializedEngineWorkflowDescriptor(id = workflowId, workflowSources = test1Sources)

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, test1Sources, Nil, Nil, localBackend)
        executionAndAuxes <- dataAccess.getWorkflowExecutionAndAuxByState(Seq(WorkflowSubmitted)) map { results =>
          results shouldNot be(empty)

          val workflowResultOption = results.find(_.execution.workflowExecutionUuid == workflowId.toString)
          workflowResultOption shouldNot be(empty)
          val workflowResult = workflowResultOption.get
          workflowResult.aux.wdlSource.toRawString should be("workflow test1 {}")
          workflowResult.aux.jsonInputs.toRawString should be("{}")
        }
      } yield ()).futureValue
    }

    it should "store and retrieve an empty String as a WdlValue" taggedAs DbmsTest in {
      val workflowId = WorkflowId.randomId()
      val workflowInfo = createMaterializedEngineWorkflowDescriptor(id = workflowId, workflowSources = test1Sources)
      val call = mock[Call]
      call.fullyQualifiedName returns "wf.a"
      val dbKey = ExecutionDatabaseKey(call.fullyQualifiedName, None, 1)
      val outputs: JobOutputs = Map("wf.a.empty" -> JobOutput(WdlString(""), None))

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, test1Sources, Nil, Nil, localBackend)
        _ <- dataAccess.setOutputs(workflowId, dbKey, outputs, Seq.empty)
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

    def publishMetadataEvents(baseKey: MetadataKey, keyValues: Array[(String, String)]): Future[Unit] = {
      val events = keyValues map { case (k, v) =>
        MetadataEvent(baseKey.copy(key = k), MetadataValue(v))
      }
      dataAccess.addMetadataEvents(events)
    }

    def baseWorkflowMetadata(name: String): Future[WorkflowId] = {
      val workflowId = WorkflowId.randomId()
      val workflowKey = MetadataKey(workflowId, jobKey = None, key = null)
      def keyAndValue(name: String) = Array(
        (WorkflowMetadataKeys.StartTime, OffsetDateTime.now.toString),
        (WorkflowMetadataKeys.Status, WorkflowSubmitted.toString),
        (WorkflowMetadataKeys.Name, name)
      )
      publishMetadataEvents(workflowKey, keyAndValue(name)).map(_ => workflowId)
    }

    it should "return pagination metadata only when page and pagesize query params are specified" taggedAs DbmsTest in {
      (for {
        workflow1Id <- baseWorkflowMetadata(Workflow1Name)
        //get metadata when page and pagesize are specified
        _ <- dataAccess.queryWorkflowSummaries(WorkflowQueryParameters(Seq(
          WorkflowQueryKey.Page.name -> "1", WorkflowQueryKey.PageSize.name -> "50"))) map { case (response, meta) =>
          meta match {
            case Some(metadata) =>
            case None => fail("Should have metadata when page and pagesize are specified.")
          }
        }
        //don't get metadata when page and pagesize are not specified
        _ <- dataAccess.queryWorkflowSummaries(
          WorkflowQueryParameters(Seq())) map { case(response, meta) =>
          meta match {
            case Some(metadata) => fail("Should not have metadata when page and pagesize are not specified")
            case None =>
          }
        }
      } yield()).futureValue
    }

    it should "create and query a workflow" taggedAs DbmsTest in {

      val randomIds = Seq.fill(10)(WorkflowId.randomId().toString)

      def succeededWorkflowMetadata(id: WorkflowId): Future[Unit] = {
        val workflowKey = MetadataKey(id, jobKey = None, key = null)
        val keyAndValue = Array(
          (WorkflowMetadataKeys.Status, WorkflowRunning.toString),
          (WorkflowMetadataKeys.Status, WorkflowSucceeded.toString),
          (WorkflowMetadataKeys.EndTime, OffsetDateTime.now.toString))

        publishMetadataEvents(workflowKey, keyAndValue)
      }

      (for {
        workflow1Id <- baseWorkflowMetadata(Workflow1Name)
        _ <- succeededWorkflowMetadata(workflow1Id)
        // Put a bit of space between the two workflows
        _ = Thread.sleep(50)
        workflow2Id <- baseWorkflowMetadata(Workflow2Name)

        // refresh the metadata
        _ <- dataAccess.refreshWorkflowMetadataSummaries(0L, None) map { max =>
          max should be > 0L
        }

        // Query with no filters
        (workflowQueryResult, workflowQueryResult2) <-
        dataAccess.queryWorkflowSummaries(WorkflowQueryParameters(Seq.empty)) map { case (response, meta) =>
          val result = response.results find { r => r.name.contains(Workflow1Name) && r.end.isDefined } getOrElse
            fail(s"$Workflow1Name with an end not found in ${response.results}")
          val result2 = response.results find {
            _.name.contains(Workflow2Name)
          } getOrElse fail(s"$Workflow2Name not found in ${response.results}")
          (result, result2)
        }
        // Filter by name
        _ <- dataAccess.queryWorkflowSummaries(WorkflowQueryParameters(Seq(WorkflowQueryKey.Name.name -> Workflow1Name))) map { case (response, meta) =>
          val resultsByName = response.results groupBy {
            _.name
          }
          resultsByName.keys.toSet.flatten should equal(Set(Workflow1Name))
        }
        // Filter by multiple names
        _ <- dataAccess.queryWorkflowSummaries(WorkflowQueryParameters(Seq(WorkflowQueryKey.Name.name -> Workflow1Name, WorkflowQueryKey.Name.name -> Workflow2Name))) map { case (response, meta) =>
          val resultsByName = response.results groupBy {
            _.name
          }
          resultsByName.keys.toSet.flatten should equal(Set(Workflow1Name, Workflow2Name))
        }
        // Filter by workflow id
        _ <- dataAccess.queryWorkflowSummaries(WorkflowQueryParameters(
          Seq(WorkflowQueryKey.Id.name -> workflow1Id.toString))) map { case (response, meta) =>
          val resultsById = response.results groupBy {
            _.name
          }
          resultsById.keys.toSet.flatten should equal(Set(Workflow1Name))
        }
        // Filter by multiple workflow ids
        _ <- dataAccess.queryWorkflowSummaries(WorkflowQueryParameters(
          Seq(workflow1Id, workflow2Id).map(id => WorkflowQueryKey.Id.name -> id.toString))) map { case (response, meta) =>
          val resultsById = response.results groupBy {
            _.name
          }
          resultsById.keys.toSet.flatten should equal(Set(Workflow1Name, Workflow2Name))
        }
        // Filter by workflow id within random Ids
        _ <- dataAccess.queryWorkflowSummaries(WorkflowQueryParameters(
          (randomIds :+ workflow1Id).map(id => WorkflowQueryKey.Id.name -> id.toString))) map { case (response, meta) =>
          val resultsById = response.results groupBy {
            _.name
          }
          resultsById.keys.toSet.flatten should equal(Set(Workflow1Name))
        }
        // Filter by status
        _ <- dataAccess.queryWorkflowSummaries(WorkflowQueryParameters(Seq(WorkflowQueryKey.Status.name -> "Submitted"))) map { case (response, meta) =>
          val resultsByStatus = response.results groupBy (_.status)
          resultsByStatus.keys.toSet.flatten should equal(Set("Submitted"))
        }
        // Filter by multiple statuses
        _ <- dataAccess.queryWorkflowSummaries(WorkflowQueryParameters(Seq(WorkflowQueryKey.Status.name -> "Submitted", WorkflowQueryKey.Status.name -> "Succeeded"))) map { case (response, meta) =>
          val resultsByStatus = response.results groupBy (_.status)
          resultsByStatus.keys.toSet.flatten should equal(Set("Submitted", "Succeeded"))
        }
        // Filter by start date
        _ <- dataAccess.queryWorkflowSummaries(WorkflowQueryParameters(Seq(
          WorkflowQueryKey.StartDate.name -> workflowQueryResult2.start.get.toString))) map { case (response, meta) =>
          response.results partition { r => r.start.isDefined && r.start.get.compareTo(workflowQueryResult.start.get) >= 0 } match {
            case (y, n) if y.nonEmpty && n.isEmpty => // good
            case (y, n) => fail(s"Found ${y.size} later workflows and ${n.size} earlier")
          }
        }
        // Filter by end date
        _ <- dataAccess.queryWorkflowSummaries(WorkflowQueryParameters(Seq(
          WorkflowQueryKey.EndDate.name -> workflowQueryResult.end.get.toString))) map { case (response, meta) =>
          response.results partition { r => r.end.isDefined && r.end.get.compareTo(workflowQueryResult.end.get) <= 0 } match {
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

    it should "support call caching configuration for specified calls in a regular workflow" taggedAs DbmsTest in {
      val sources = SampleWdl.ThreeStep.asWorkflowSources()
      val id = WorkflowId.randomId()
      val workflowInfo = createMaterializedEngineWorkflowDescriptor(id, sources)

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, sources, Nil, workflowInfo.namespace.workflow.calls, localBackend)
        // Unknown workflow
        _ <- assertCallCachingFailure(WorkflowId.randomId(), callName = Option("three_step.ps"), "Workflow not found")
        _ <- dataAccess.setTerminalStatus(workflowInfo.id, ExecutionDatabaseKey("three_step.ps", None, 1), ExecutionStatus.Done, None, None, None)
        executions <- dataAccess.getExecutions(workflowInfo.id.toString)
        _ = executions should have size 3
        _ = executions foreach { _.allowsResultReuse shouldBe true }
        params <- CallCachingParameters.from(workflowInfo.id, Option("three_step.ps"), AllowFalse, dataAccess)
        _ = params.workflowId shouldBe workflowInfo.id
        _ = params.callKey shouldEqual Option(ExecutionDatabaseKey("three_step.ps", None, 1))
        _ <- dataAccess.updateCallCaching(params)
        executions <- dataAccess.getExecutions(workflowInfo.id.toString)
        (allowing, disallowing) = executions partition { _.allowsResultReuse }
        _ = allowing should have size 2
        _ = disallowing should have size 1
        _ = disallowing.seq.head.callFqn should be("three_step.ps")
      } yield ()).futureValue
    }

    it should "support call caching configuration for specified calls in a scattered workflow" taggedAs DbmsTest in {
      val sources = SampleWdl.SimpleScatterWdl.asWorkflowSources()
      val id = WorkflowId.randomId()
      val workflowInfo = createMaterializedEngineWorkflowDescriptor(id, sources)

      (for {
      // The `inside_scatter` is a to-be-exploded placeholder, but it will conflict with the collector that the
      // scatter explodes below so filter that out.
        _ <- dataAccess.createWorkflow(workflowInfo, sources, Nil, workflowInfo.namespace.workflow.calls.filterNot(_.unqualifiedName == "inside_scatter"), localBackend)
        scatter = workflowInfo.namespace.workflow.scatters.head
        scatterKey = ScatterKey(scatter, None)
        newEntries = scatterKey.populate(5)
        _ <- dataAccess.insertCalls(workflowInfo.id, newEntries.keys, localBackend)
        executions <- dataAccess.getExecutions(workflowInfo.id.toString)
        _ = executions foreach { _.allowsResultReuse shouldBe true }

        // Calls outside the scatter should work the same as an unscattered workflow.
        outsideParams <- CallCachingParameters.from(workflowInfo.id, Option("scatter0.outside_scatter"), AllowFalse, dataAccess)
        _ <- dataAccess.updateCallCaching(outsideParams)
        executions <- dataAccess.getExecutions(workflowInfo.id.toString)
        (allowing, disallowing) = executions partition { _.allowsResultReuse }
        _ = allowing should have size (executions.size - 1)
        _ = disallowing should have size 1
        _ = disallowing.seq.head.callFqn should be ("scatter0.outside_scatter")

        // Support unindexed scattered call targets to update all shards.
        _ = executions filter { _.callFqn == "scatter0.inside_scatter" } foreach { _.allowsResultReuse shouldBe true }
        unindexedCallParams <- CallCachingParameters.from(workflowInfo.id, Option("scatter0.inside_scatter"), AllowFalse, dataAccess)
        _ <- dataAccess.updateCallCaching(unindexedCallParams)
        executions <- dataAccess.getExecutions(workflowInfo.id.toString)
        _ = executions filter { _.callFqn == "scatter0.inside_scatter" } foreach { _.allowsResultReuse shouldBe false }

        // Support indexed shards as well.
        insideParams <- CallCachingParameters.from(workflowInfo.id, Option("scatter0.inside_scatter.3"), AllowTrue, dataAccess)
        _ <- dataAccess.updateCallCaching(insideParams)
        executions <- dataAccess.getExecutions(workflowInfo.id.toString)
        inside = executions filter { e => e.callFqn == "scatter0.inside_scatter" }
        (allowing, disallowing) = inside partition { _.allowsResultReuse }
        _ = allowing should have size 1
        _ = disallowing should have size (inside.size - 1)
      } yield ()).futureValue
    }

    it should "support call caching configuration for all calls in a regular workflow" taggedAs DbmsTest in {
      val sources = SampleWdl.ThreeStep.asWorkflowSources()
      val id = WorkflowId.randomId()
      val workflowInfo = createMaterializedEngineWorkflowDescriptor(id, sources)
      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, sources, Nil, workflowInfo.namespace.workflow.calls, localBackend)
        // Unknown workflow
        _ <- assertCallCachingFailure(WorkflowId.randomId(), callName = None, "Workflow not found")
        params <- CallCachingParameters.from(workflowInfo.id, None, AllowFalse, dataAccess)
        executions <- dataAccess.getExecutions(workflowInfo.id.toString)
        _ = executions should have size 3
        _ = executions foreach { _.allowsResultReuse shouldBe true }
        _ <- dataAccess.updateCallCaching(params)
        executions <- dataAccess.getExecutions(workflowInfo.id.toString)
        _ = executions foreach { _.allowsResultReuse shouldBe false }
      } yield ()).futureValue
    }

    it should "support call caching configuration for all calls in a scattered workflow" taggedAs DbmsTest in {
      val sources = SampleWdl.SimpleScatterWdl.asWorkflowSources()
      val id = WorkflowId.randomId()
      val workflowInfo = createMaterializedEngineWorkflowDescriptor(id, sources)
      (for {
      // The `inside_scatter` is a to-be-exploded placeholder, but it will conflict with the collector that the
      // scatter explodes below so filter that out.
        _ <- dataAccess.createWorkflow(workflowInfo, sources, Nil, workflowInfo.namespace.workflow.calls.filterNot(_.unqualifiedName == "inside_scatter"), localBackend)
        scatter = workflowInfo.namespace.workflow.scatters.head
        scatterKey = ScatterKey(scatter, None)
        newEntries = scatterKey.populate(5)
        _ <- dataAccess.insertCalls(workflowInfo.id, newEntries.keys, localBackend)

        scatterParams <- CallCachingParameters.from(workflowInfo.id, None, AllowFalse, dataAccess)
        executions <- dataAccess.getExecutions(workflowInfo.id.toString)
        _ = executions foreach { _.allowsResultReuse shouldBe true }
        _ <- dataAccess.updateCallCaching(scatterParams)
        executions <- dataAccess.getExecutions(workflowInfo.id.toString)
        _ = executions foreach { _.allowsResultReuse shouldBe false }
      } yield ()).futureValue
    }


    it should "query a single execution status" taggedAs DbmsTest in {
      val workflowId = WorkflowId.randomId()
      val workflowInfo = createMaterializedEngineWorkflowDescriptor(id = workflowId, workflowSources = test1Sources)
      val task = Task.empty
      val callFqn = "fully.qualified.name"
      val call = new Call(None, callFqn, task, Set.empty[FullyQualifiedName], Map.empty, None)

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, test1Sources, Nil, Seq(call), localBackend)
        _ <- dataAccess.updateWorkflowState(workflowId, WorkflowRunning)
        _ <- dataAccess.getExecutionStatus(workflowId, ExecutionDatabaseKey(callFqn, None, 1)) map { status =>
          status.get.executionStatus shouldBe ExecutionStatus.NotStarted
          status.get.returnCode shouldBe None
        }
      } yield ()).futureValue
    }

    it should "create and retrieve 3step.wdl with a 10,000 char pattern" taggedAs DbmsTest in {
      val workflowId = WorkflowId.randomId()
      val sampleWdl = SampleWdl.ThreeStepLargeJson
      val sources = sampleWdl.asWorkflowSources()
      val workflowInfo = createMaterializedEngineWorkflowDescriptor(id = workflowId, workflowSources = sources)

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, sources, Nil, Nil, localBackend)
        executionAndAuxes <- dataAccess.getWorkflowExecutionAndAuxByState(Seq(WorkflowSubmitted)) map { results =>
          results shouldNot be(empty)

          val workflowResultOption = results.find(_.execution.workflowExecutionUuid == workflowId.toString)
          workflowResultOption shouldNot be(empty)
          val workflowResult = workflowResultOption.get
          workflowResult.aux.wdlSource.toRawString should be(sampleWdl.wdlSource())
          workflowResult.aux.jsonInputs.toRawString should be(sampleWdl.wdlJson)
        }
      } yield ()).futureValue
    }

    it should "fail when saving a workflow twice" taggedAs DbmsTest in {
      val workflowId = WorkflowId.randomId()
      val workflowInfo = createMaterializedEngineWorkflowDescriptor(id = workflowId, workflowSources = test1Sources)

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, test1Sources, Nil, Nil, localBackend)
        _ <- dataAccess.createWorkflow(workflowInfo, test1Sources, Nil, Nil, localBackend)
      } yield ()).failed.futureValue should be(a[SQLException])
    }

    it should "fail when updating a non-existent workflow state" taggedAs DbmsTest in {
      val workflowId = WorkflowId.randomId()

      (for {
        _ <- dataAccess.updateWorkflowState(workflowId, WorkflowRunning)
      } yield ()).failed.futureValue should be(an[RuntimeException])
    }

    it should "update and get a workflow state" taggedAs DbmsTest in {
      val workflowId = WorkflowId.randomId()
      val workflowInfo = createMaterializedEngineWorkflowDescriptor(id = workflowId, workflowSources = test1Sources)

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, test1Sources, Nil, Nil, localBackend)
        _ <- dataAccess.updateWorkflowState(workflowId, WorkflowRunning)
        _ <- dataAccess.getWorkflowExecutionAndAuxByState(Seq(WorkflowRunning)) map { results =>
          results shouldNot be(empty)

          val workflowResultOption = results.find(_.execution.workflowExecutionUuid == workflowId.toString)
          workflowResultOption shouldNot be(empty)
          val workflowResult = workflowResultOption.get
          workflowResult.aux.wdlSource.toRawString should be("workflow test1 {}")
          workflowResult.aux.jsonInputs.toRawString should be("{}")
        }
      } yield ()).futureValue
    }

    it should "get workflow state" taggedAs DbmsTest in {
      val workflowId = WorkflowId.randomId()
      val workflowInfo = createMaterializedEngineWorkflowDescriptor(id = workflowId, workflowSources = test1Sources)
      val key = new SymbolStoreKey("myScope", "myName", None, input = true)
      val entry = new SymbolStoreEntry(key, WdlStringType, Option(testWdlString), Option(testWdlStringHash))

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, test1Sources, Seq(entry), Nil, localBackend)
        _ <- dataAccess.getWorkflowExecution(workflowId.toString) map { w =>
          w.status shouldBe WorkflowSubmitted.toString
          w.endDt shouldBe None
        }
        _ <- dataAccess.updateWorkflowState(workflowId, WorkflowSucceeded)
        _ <- dataAccess.getWorkflowExecution(workflowId.toString) map { w =>
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

      it should s"get $spec" taggedAs DbmsTest in {
        val workflowId = WorkflowId.randomId()
        val workflowInfo = createMaterializedEngineWorkflowDescriptor(id = workflowId, workflowSources = test1Sources)

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
          _ <- dataAccess.createWorkflow(workflowInfo, test1Sources, Nil, Seq(call), localBackend)
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

    it should "insert a call in execution table" taggedAs DbmsTest in {
      val callFqn = "call.fully.qualified.scope"
      val workflowId = WorkflowId.randomId()
      val workflowInfo = createMaterializedEngineWorkflowDescriptor(id = workflowId, workflowSources = test1Sources)
      val task = Task.empty
      val call = new Call(None, callFqn, task, Set.empty[FullyQualifiedName], Map.empty, None)
      val callKey = BackendCallKey(call, None, 1)
      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, test1Sources, Nil, Nil, localBackend)
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

    it should "fail to get an non-existent execution status" taggedAs DbmsTest in {
      dataAccess.getExecutionStatuses(WorkflowId.randomId()).failed.futureValue should be(a[NoSuchElementException])
    }

    it should "get a symbol input" taggedAs DbmsTest in {
      val callFqn = "call.fully.qualified.scope"
      val symbolFqn = "symbol.fully.qualified.scope"
      val workflowId = WorkflowId.randomId()
      val workflowInfo = createMaterializedEngineWorkflowDescriptor(id = workflowId, workflowSources = test1Sources)
      val key = new SymbolStoreKey(callFqn, symbolFqn, None, input = true)
      val entry = new SymbolStoreEntry(key, WdlStringType, Option(testWdlString), Option(testWdlStringHash))
      val task = Task.empty
      val call = new Call(None, callFqn, task, Set.empty[FullyQualifiedName], Map.empty, None)

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, test1Sources, Seq(entry), Nil, localBackend)
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

    it should "get a symbol input that has a very long WDL value field" taggedAs DbmsTest in {
      val wdlArray = new WdlArray(WdlArrayType(WdlStringType), Seq(WdlString(Workflow1Name), WdlString("*" * 10000)))
      val callFqn = "call.fully.qualified.scope"
      val symbolFqn = "symbol.fully.qualified.scope"
      val workflowId = WorkflowId.randomId()
      val workflowDescriptor = createMaterializedEngineWorkflowDescriptor(id = workflowId, workflowSources = test1Sources)
      val key = new SymbolStoreKey(callFqn, symbolFqn, None, input = true)
      // PBE hack
      // val entry = new SymbolStoreEntry(key, WdlArrayType(WdlStringType), Option(wdlArray), workflowDescriptor.hash(wdlArray))
      val junkEmptyHash = None
      val entry = new SymbolStoreEntry(key, WdlArrayType(WdlStringType), Option(wdlArray), junkEmptyHash)
      val task = Task.empty
      val call = new Call(None, callFqn, task, Set.empty[FullyQualifiedName], Map.empty, None)

      (for {
        _ <- dataAccess.createWorkflow(workflowDescriptor, test1Sources, Seq(entry), Nil, localBackend)
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
          resultSymbol.symbolHash should be(empty)
        }
      } yield ()).futureValue
    }

    it should "fail to get inputs for a null call" taggedAs DbmsTest in {
      val workflowId: WorkflowId = WorkflowId.randomId()
      val workflowInfo = createMaterializedEngineWorkflowDescriptor(id = workflowId, workflowSources = test1Sources)

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, test1Sources, Nil, Nil, localBackend)
        _ <- dataAccess.updateWorkflowState(workflowId, WorkflowRunning)
        _ <- dataAccess.getInputs(workflowId, null)
      } yield ()).failed.futureValue should be(an[RuntimeException])
    }

    it should "set and get an output" taggedAs DbmsTest in {
      val callFqn = "call.fully.qualified.scope"
      val symbolLqn = "symbol"
      val workflowId = WorkflowId.randomId()
      val workflowInfo = createMaterializedEngineWorkflowDescriptor(id = workflowId, workflowSources = test1Sources)
      val task = Task.empty
      val call = new Call(None, callFqn, task, Set.empty[FullyQualifiedName], Map.empty, None)

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, test1Sources, Nil, Nil, localBackend)
        _ <- dataAccess.updateWorkflowState(workflowId, WorkflowRunning)
        _ <- dataAccess.setOutputs(workflowId, ExecutionDatabaseKey(call.fullyQualifiedName, None, 1),
          Map(symbolLqn -> JobOutput(testWdlString, Option(testWdlStringHash))), Seq.empty)
        _ <- dataAccess.setOutputs(workflowId, ExecutionDatabaseKey(call.fullyQualifiedName, Option(0), 1),
          Map(symbolLqn -> JobOutput(testWdlStringShard, Option(testWdlStringShardHash))), Seq.empty)
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

    it should "set and get shard statuses" taggedAs DbmsTest in {
      val callFqn1 = "call.fully.qualified.scope$s1"
      val callFqn2 = "call.fully.qualified.scope$s2"
      val workflowId = WorkflowId.randomId()
      val workflowInfo = createMaterializedEngineWorkflowDescriptor(id = workflowId, workflowSources = test1Sources)
      val task = Task.empty
      val call1 = new Call(None, callFqn1, task, Set.empty[FullyQualifiedName], Map.empty, None)
      val shardIndex1 = Option(0)
      val pid1 = Option("123")
      val call2 = new Call(None, callFqn2, task, Set.empty[FullyQualifiedName], Map.empty, None)
      val shardIndex2 = Option(1)
      val pid2 = Option("987")

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, test1Sources, Nil, Nil, localBackend)
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

    it should "set and get an output by call" taggedAs DbmsTest in {
      val callFqn = "call.fully.qualified.scope"
      val symbolLqn = "symbol"
      val workflowId = WorkflowId.randomId()
      val workflowInfo = createMaterializedEngineWorkflowDescriptor(id = workflowId, workflowSources = test1Sources)
      val task = Task.empty
      val call = new Call(None, callFqn, task, Set.empty[FullyQualifiedName], Map.empty, None)

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, test1Sources, Nil, Nil, localBackend)
        _ <- dataAccess.updateWorkflowState(workflowId, WorkflowRunning)
        _ <- dataAccess.setOutputs(workflowId, ExecutionDatabaseKey(call.fullyQualifiedName, None, 1),
          Map(symbolLqn -> JobOutput(testWdlString, Option(testWdlStringHash))), Seq.empty)
        _ <- dataAccess.setOutputs(workflowId, ExecutionDatabaseKey(call.fullyQualifiedName, Option(0), 1),
          Map(symbolLqn -> JobOutput(testWdlStringShard, Option(testWdlStringShardHash))), Seq.empty)
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

    it should "fail to get outputs for a null call" taggedAs DbmsTest in {
      val workflowId = WorkflowId.randomId()
      val workflowInfo = createMaterializedEngineWorkflowDescriptor(id = workflowId, workflowSources = test1Sources)

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, test1Sources, Nil, Nil, localBackend)
        _ <- dataAccess.updateWorkflowState(workflowId, WorkflowRunning)
        _ <- dataAccess.getOutputs(workflowId, null)
      } yield ()).failed.futureValue should be(an[RuntimeException])
    }

    it should "fail to create workflow for an unknown backend" taggedAs DbmsTest in {
      val callFqn = "call.fully.qualified.scope"
      val workflowId = WorkflowId.randomId()
      val workflowInfo = createMaterializedEngineWorkflowDescriptor(id = workflowId, workflowSources = test1Sources)
      val task = Task.empty
      val call = new Call(None, callFqn, task, Set.empty[FullyQualifiedName], Map.empty, None)

      dataAccess.createWorkflow(workflowInfo, test1Sources, Nil, Seq(call),
        UnknownOldStyleBackend$).failed.futureValue should be(an[Exception])
    }

    it should "fail to create workflow for a null backend" taggedAs DbmsTest in {
      val callFqn = "call.fully.qualified.scope"
      val workflowId = WorkflowId.randomId()
      val workflowInfo = createMaterializedEngineWorkflowDescriptor(id = workflowId, workflowSources = test1Sources)
      val task = Task.empty
      val call = new Call(None, callFqn, task, Set.empty[FullyQualifiedName], Map.empty, None)

      dataAccess.createWorkflow(workflowInfo, test1Sources, Nil, Seq(call),
        null).failed.futureValue should be(a[NullPointerException])
    }

    it should "set and get the same symbol with IO as input then output" taggedAs DbmsTest in {
      val callFqn = "call.fully.qualified.scope"
      val symbolLqn = "symbol"
      val symbolFqn = callFqn + "." + symbolLqn
      val workflowId = WorkflowId.randomId()
      val workflowInfo = createMaterializedEngineWorkflowDescriptor(id = workflowId, workflowSources = test1Sources)
      val key = new SymbolStoreKey(callFqn, symbolFqn, None, input = true)
      val entry = new SymbolStoreEntry(key, WdlStringType, Option(testWdlString), Option(testWdlStringHash))
      val task = Task.empty
      val call = new Call(None, callFqn, task, Set.empty[FullyQualifiedName], Map.empty, None)

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, test1Sources, Seq(entry), Nil, localBackend)
        _ <- dataAccess.updateWorkflowState(workflowId, WorkflowRunning)
        _ <- dataAccess.setOutputs(workflowId, ExecutionDatabaseKey(call.fullyQualifiedName, None, 1),
          Map(symbolLqn -> JobOutput(testWdlString, Option(testWdlStringHash))), Seq.empty)
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

    it should "fail when setting an existing symbol output" taggedAs DbmsTest in {
      val callFqn = "call.fully.qualified.scope"
      val symbolLqn = "symbol"
      val symbolFqn = callFqn + "." + symbolLqn
      val workflowId = WorkflowId.randomId()
      val workflowInfo = createMaterializedEngineWorkflowDescriptor(id = workflowId, workflowSources = test1Sources)
      val task = Task.empty
      val call = new Call(None, callFqn, task, Set.empty[FullyQualifiedName], Map.empty, None)

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, test1Sources, Seq(), Nil, localBackend)
        _ <- dataAccess.updateWorkflowState(workflowId, WorkflowRunning)
        _ <- dataAccess.setOutputs(workflowId, ExecutionDatabaseKey(call.fullyQualifiedName, None, 1),
          Map(symbolFqn -> JobOutput(testWdlString, Option(testWdlStringHash))), Seq.empty)
        // Second attempt should fail
        _ <- dataAccess.setOutputs(workflowId, ExecutionDatabaseKey(call.fullyQualifiedName, None, 1),
          Map(symbolFqn -> JobOutput(testWdlString, Option(testWdlStringHash))), Seq.empty)
      } yield ()).failed.futureValue should be(a[SQLException])
    }

    it should "set and get a backend info" taggedAs DbmsTest in {
      val workflowId = WorkflowId.randomId()
      val workflowInfo = createMaterializedEngineWorkflowDescriptor(id = workflowId, workflowSources = test1Sources)
      val task = Task.empty
      val call = new Call(None, "fully.qualified.name", task, Set.empty[FullyQualifiedName], Map.empty, None)
      val pid = Option("123")

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo, test1Sources, Nil, Seq(call), localBackend)
        _ <- dataAccess.updateExecutionInfo(workflowId, BackendCallKey(call, None, 1), InfoKeys.Pid, pid)
        _ <- dataAccess.getExecutionInfos(workflowId, call, 1) map { info =>
          info should have size 1
          info.head.key should equal(InfoKeys.Pid)
          info.head.value should equal(pid)
        }
      } yield ()).futureValue
    }

    // Queries use `.head` a lot. There was a bug that pulled the backend info by fqn, but for any workflow.
    it should "set and get a backend info for same call on two workflows" taggedAs DbmsTest in {
      val workflowId1 = WorkflowId.randomId()
      val workflowId2 = WorkflowId.randomId()
      val workflowInfo1 = createMaterializedEngineWorkflowDescriptor(id = workflowId1, workflowSources = test1Sources)
      val workflowInfo2 = createMaterializedEngineWorkflowDescriptor(id = workflowId2, workflowSources = test1Sources)
      val task = Task.empty
      val call = new Call(None, "fully.qualified.name", task, Set.empty[FullyQualifiedName], Map.empty, None)
      val pid1 = Option("123")
      val pid2 = Option("321")

      (for {
        _ <- dataAccess.createWorkflow(workflowInfo1, test1Sources, Nil, Seq(call), localBackend)
        _ <- dataAccess.createWorkflow(workflowInfo2, test1Sources, Nil, Seq(call), localBackend)
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

    it should "get execution infos by key" taggedAs DbmsTest in {
      val workflowId1 = WorkflowId.randomId()
      val sources = WorkflowSourceFiles("task a {command {}} workflow test1 {call a}", "{}", "{}")
      val descriptor = createMaterializedEngineWorkflowDescriptor(
        id = workflowId1,
        workflowSources = sources
      )
      val call = descriptor.namespace.workflow.calls.find(_.unqualifiedName == "a").get
      val infos = Map(
        "key0" -> Option("foo"),
        "key1" -> None
      )

      (for {
        _ <- dataAccess.createWorkflow(descriptor, sources, Nil, Seq(call), localBackend)
        _ <- dataAccess.upsertExecutionInfo(workflowId1, BackendCallKey(call, None, 1), infos, actorSystem)
        _ <- dataAccess.getExecutionInfoByKey(workflowId1, call, 1, "key0") map { info =>
          info shouldEqual Some(Some("foo"))
        }
        _ <- dataAccess.getExecutionInfoByKey(workflowId1, call, 1, "key1") map { info =>
          info shouldEqual Some(None)
        }
        _ <- dataAccess.getExecutionInfoByKey(workflowId1, call, 1, "bad") map { info =>
          info shouldEqual None
        }
      } yield ()).futureValue
    }

    it should "update call start and end dates appropriately" taggedAs DbmsTest in {
      val callFqn = "call.fully.qualified.scope"
      val workflowId = WorkflowId.randomId()
      val workflowInfo = createMaterializedEngineWorkflowDescriptor(id = workflowId, workflowSources = test1Sources)
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
        _ <- dataAccess.createWorkflow(workflowInfo, test1Sources, Nil, Seq(call), localBackend)
        callKey = ExecutionDatabaseKey(callFqn, None, 1)

        _ <- dataAccess.updateStatus(workflowId, List(callKey), ExecutionStatus.NotStarted)
        _ <- dataAccess.getExecutions(workflowId.toString) map assertDates(startDefined = false, endDefined = false)

        _ <- dataAccess.setStartingStatus(workflowId, List(callKey))
        _ <- dataAccess.getExecutions(workflowId.toString) map assertDates(startDefined = true, endDefined = false)

        _ <- dataAccess.setTerminalStatus(workflowId, callKey, ExecutionStatus.Done, None, None, None)
        _ <- dataAccess.getExecutions(workflowId.toString) map assertDates(startDefined = true, endDefined = true)
      } yield()).futureValue
    }

    it should "not deadlock with upserts" taggedAs DbmsTest in {
      val workflowId = WorkflowId.randomId()
      val sources = WorkflowSourceFiles(
        wdlSource="""
          |task a {command{}}
          |workflow w {
          |  call a
          |  call a as b
          |  call a as c
          |}
        """.stripMargin,
        inputsJson="{}",
        workflowOptionsJson="{}"
      )
      val descriptor = createMaterializedEngineWorkflowDescriptor(id = workflowId, workflowSources = sources)
      val key1 = ExecutionDatabaseKey("w.a", None, 1)
      val key2 = ExecutionDatabaseKey("w.b", None, 1)
      val key3 = ExecutionDatabaseKey("w.c", None, 1)
      val attributes = Map(
        "a" -> WdlString("foo"),
        "b" -> WdlInteger(1),
        "c" -> WdlFloat(2.2),
        "d" -> WdlString("foo"),
        "e" -> WdlInteger(1),
        "f" -> WdlFloat(2.2)
      )
      (for {
        _ <- dataAccess.createWorkflow(descriptor, sources, Nil, descriptor.namespace.workflow.calls, localBackend)
        executions <- dataAccess.getExecutions(descriptor.id.toString)
        _ = executions should have size 3
        _ <- Future.sequence(Seq(
          dataAccess.upsertRuntimeAttributes(workflowId, key1, attributes, actorSystem),
          dataAccess.upsertRuntimeAttributes(workflowId, key2, attributes, actorSystem),
          dataAccess.upsertRuntimeAttributes(workflowId, key3, attributes, actorSystem)
        ))
      } yield ()).futureValue
    }

    it should "close the database" taggedAs DbmsTest in {
      dataAccess.close()
    }
  }
}
