package cromwell.engine.db.slick

import java.nio.file.{FileSystem, FileSystems}
import java.sql.SQLException
import java.time.OffsetDateTime

import better.files._
import com.google.api.client.util.ExponentialBackOff
import com.typesafe.config.ConfigFactory
import cromwell.CromwellSpec.DbmsTest
import cromwell.CromwellTestkitSpec.TestWorkflowManagerSystem
import cromwell.core._
import cromwell.database.SqlConverters._
import cromwell.database.obj.{Execution, WorkflowMetadataKeys}
import cromwell.database.slick.SlickDatabase
import cromwell.engine._
import cromwell.engine.backend._
import cromwell.engine.backend.io._
import cromwell.engine.db.DataAccess
import cromwell.services.{MetadataEvent, MetadataKey, MetadataValue}
import cromwell.util.SampleWdl
import cromwell.webservice.{WorkflowQueryKey, WorkflowQueryParameters}
import cromwell.{CromwellTestkitSpec, webservice}
import org.scalactic.StringNormalizations._
import org.scalatest.concurrent.PatienceConfiguration.Timeout
import org.scalatest.concurrent.ScalaFutures
import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.time.SpanSugar._
import org.scalatest.time.{Millis, Seconds, Span}
import org.scalatest.{BeforeAndAfterAll, FlatSpec, Matchers}
import org.specs2.mock.Mockito
import wdl4s._
import wdl4s.types.{WdlArrayType, WdlStringType}
import wdl4s.values._

import scala.concurrent.{ExecutionContext, Future}
import scala.util.Try

object SlickDataAccessSpec {
  val AllowFalse = Seq(webservice.QueryParameter("allow", "false"))
  val AllowTrue = Seq(webservice.QueryParameter("allow", "true"))

  val Workflow1Name = "test1"
  val Workflow2Name = "test2"

//TODO: Why is this still here?
//  implicit class EnhancedEngineWorkflowDescriptor(val workflowDescriptor: EngineWorkflowDescriptor) extends AnyVal {
//    def fileHasher: FileHasher = { wdlFile: WdlFile =>
//        SymbolHash(wdlFile.value.toPath(FileSystems.getDefault).hash)
//    }
//  }
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

  //lazy val localBackend = OldStyleLocalBackend(CromwellTestkitSpec.DefaultLocalBackendConfigEntry, workflowManagerSystem.actorSystem)

  val test1Sources = WorkflowSourceFiles("workflow test1 {}", "{}", "{}")
  val test2Sources = WorkflowSourceFiles("workflow test2 {}", "{}", "{}")

  it should "not deadlock" ignore {
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
//    def assertTerminalExecution(id: WorkflowId, expectedCount: Int, expectedTerminal: ExpectedTerminal): Future[Unit] = {
//      for {
//        _ <- dataAccess.getExecutions(id.toString) map { executions =>
//          executions should have size expectedCount
//          executions foreach { execution =>
//            execution.endDt.isDefined shouldBe expectedTerminal(execution.callFqn, execution.index)
//          }
//        }
//      } yield ()
//    }

    val id = WorkflowId.randomId()

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

    it should "create and retrieve the workflow for just reading" taggedAs DbmsTest ignore {
//      val workflowId = WorkflowId.randomId()
//      val workflowInfo = createMaterializedEngineWorkflowDescriptor(id = workflowId, workflowSources = test1Sources)
//
//      (for {
//        _ <- dataAccess.createWorkflow(workflowInfo, test1Sources, Nil, Nil, localBackend)
//        executionAndAuxes <- dataAccess.getWorkflowExecutionAndAuxByState(Seq(WorkflowSubmitted)) map { results =>
//          results shouldNot be(empty)
//
//          val workflowResultOption = results.find(_.execution.workflowExecutionUuid == workflowId.toString)
//          workflowResultOption shouldNot be(empty)
//          val workflowResult = workflowResultOption.get
//          workflowResult.aux.wdlSource.toRawString should be("workflow test1 {}")
//          workflowResult.aux.jsonInputs.toRawString should be("{}")
//        }
//      } yield ()).futureValue
    }

    it should "store and retrieve an empty String as a WdlValue" taggedAs DbmsTest ignore {
//      val workflowId = WorkflowId.randomId()
//      val workflowInfo = createMaterializedEngineWorkflowDescriptor(id = workflowId, workflowSources = test1Sources)
//      val call = mock[Call]
//      call.fullyQualifiedName returns "wf.a"
//      val dbKey = ExecutionDatabaseKey(call.fullyQualifiedName, None, 1)
//      val outputs: JobOutputs = Map("wf.a.empty" -> JobOutput(WdlString(""), None))
//
//      (for {
//        _ <- dataAccess.createWorkflow(workflowInfo, test1Sources, Nil, Nil, localBackend)
//        _ <- dataAccess.setOutputs(workflowId, dbKey, outputs, Seq.empty)
//        _ <- dataAccess.getOutputs(workflowId, dbKey) map { results =>
//          results shouldNot be(empty)
//
//          val workflowOutput = results.headOption
//          workflowOutput shouldNot be(empty)
//          val workflowResult = workflowOutput.get
//          workflowResult.wdlValue shouldNot be(empty)
//          val wdlValue = workflowResult.wdlValue
//          wdlValue.get.valueString shouldBe ""
//        }
//      } yield ()).futureValue
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

    it should "close the database" taggedAs DbmsTest in {
      dataAccess.close()
    }
  }
}
