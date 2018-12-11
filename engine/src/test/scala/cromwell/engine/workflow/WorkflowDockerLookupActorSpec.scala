package cromwell.engine.workflow

import akka.actor.{ActorRef, Props}
import akka.testkit.{ImplicitSender, TestActorRef, TestProbe}
import com.typesafe.config.ConfigFactory
import common.util.Backoff
import cromwell.core.actor.StreamIntegration.BackPressure
import cromwell.core.retry.SimpleExponentialBackoff
import cromwell.core.{TestKitSuite, WorkflowId}
import cromwell.database.slick.EngineSlickDatabase
import cromwell.database.sql.tables.DockerHashStoreEntry
import cromwell.docker.DockerInfoActor.{DockerInfoFailedResponse, DockerInfoSuccessResponse, DockerInformation}
import cromwell.docker.{DockerHashResult, DockerImageIdentifier, DockerImageIdentifierWithoutHash, DockerInfoRequest}
import cromwell.engine.workflow.WorkflowDockerLookupActor.{DockerHashActorTimeout, WorkflowDockerLookupFailure, WorkflowDockerTerminalFailure}
import cromwell.engine.workflow.WorkflowDockerLookupActorSpec._
import cromwell.engine.workflow.workflowstore.{StartableState, Submitted}
import cromwell.services.EngineServicesStore
import cromwell.services.ServicesStore._
import org.scalatest.{BeforeAndAfter, FlatSpecLike, Matchers}
import org.specs2.mock.Mockito

import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, Future}
import scala.language.postfixOps
import scala.util.control.NoStackTrace


class WorkflowDockerLookupActorSpec extends TestKitSuite("WorkflowDockerLookupActorSpecSystem") with FlatSpecLike with Matchers with ImplicitSender with BeforeAndAfter with Mockito {
  var workflowId: WorkflowId = _
  var dockerHashingActor: TestProbe = _
  var numReads: Int = _
  var numWrites: Int = _

  before {
    workflowId = WorkflowId.randomId()
    dockerHashingActor = TestProbe()
    numReads = 0
    numWrites = 0
  }

  it should "wait and resubmit the docker request when it gets a backpressure message" in {
    val backoff = SimpleExponentialBackoff(2.seconds, 10.minutes, 2D)

    val lookupActor = TestActorRef(Props(new TestWorkflowDockerLookupActor(workflowId, dockerHashingActor.ref, Submitted, backoff)), self)
    lookupActor ! LatestRequest

    dockerHashingActor.expectMsg(LatestRequest)
    dockerHashingActor.reply(BackPressure(LatestRequest))
    // Give a couple of seconds of margin to account for test latency etc...
    dockerHashingActor.expectMsg(2.seconds.+(5 seconds), LatestRequest)
  }

  it should "not look up the same tag again after a successful lookup" in {
    val db = dbWithWrite {
      numWrites = numWrites + 1
      Future.successful(())
    }

    val lookupActor = TestActorRef(WorkflowDockerLookupActor.props(workflowId, dockerHashingActor.ref, isRestart = false, db))
    lookupActor ! LatestRequest

    // The WorkflowDockerLookupActor should not have the hash for this tag yet and will need to query the dockerHashingActor.
    dockerHashingActor.expectMsg(LatestRequest)
    dockerHashingActor.reply(LatestSuccessResponse)
    // The WorkflowDockerLookupActor should forward the success message to this actor.
    expectMsg(LatestSuccessResponse)
    numWrites should equal(1)

    // Now the WorkflowDockerLookupActor should now have this hash in its mappings and should not query the dockerHashingActor again.
    lookupActor ! LatestRequest
    dockerHashingActor.expectNoMessage()
    // The WorkflowDockerLookupActor should forward the success message to this actor.
    expectMsg(LatestSuccessResponse)
    numWrites should equal(1)
  }

  it should "soldier on after docker hashing actor timeouts" in {
    val lookupActor = TestActorRef(WorkflowDockerLookupActor.props(workflowId, dockerHashingActor.ref, isRestart = false))

    lookupActor ! LatestRequest
    lookupActor ! OlderRequest

    val timeout = DockerHashActorTimeout(LatestRequest)

    // The WorkflowDockerLookupActor should not have the hash for this tag yet and will need to query the dockerHashingActor.
    dockerHashingActor.expectMsg(LatestRequest)
    dockerHashingActor.expectMsg(OlderRequest)
    dockerHashingActor.reply(timeout)
    // Send a response for the older request after sending the timeout.  This should cause a mapping to be entered in the
    // WorkflowDockerLookupActor for the older request, which will keep the WorkflowDockerLookupActor from querying the
    // DockerHashActor for this hash again.
    dockerHashingActor.reply(OlderSuccessResponse)

    val results = receiveN(2, 2 seconds).toSet
    val failedRequests = results collect {
      case f: WorkflowDockerLookupFailure if f.request == LatestRequest => f.request
    }

    failedRequests should equal(Set(LatestRequest))

    // Try again.  The hashing actor should receive the latest message and this time won't time out.
    lookupActor ! LatestRequest
    lookupActor ! OlderRequest
    dockerHashingActor.expectMsg(LatestRequest)
    dockerHashingActor.reply(LatestSuccessResponse)

    val responses = receiveN(2, 2 seconds).toSet
    val hashResponses = responses collect { case msg: DockerInfoSuccessResponse => msg }
    // Success after transient timeout failures:
    hashResponses should equal(Set(LatestSuccessResponse, OlderSuccessResponse))
  }

  it should "respond appropriately to docker hash lookup failures" in {
    val lookupActor = TestActorRef(WorkflowDockerLookupActor.props(workflowId, dockerHashingActor.ref, isRestart = false))
    lookupActor ! LatestRequest
    lookupActor ! OlderRequest

    // The WorkflowDockerLookupActor should not have the hash for this tag yet and will need to query the dockerHashingActor.
    dockerHashingActor.expectMsg(LatestRequest)
    dockerHashingActor.expectMsg(OlderRequest)
    val olderFailedResponse = DockerInfoFailedResponse(new RuntimeException("Lookup failed"), OlderRequest)

    dockerHashingActor.reply(LatestSuccessResponse)
    dockerHashingActor.reply(olderFailedResponse)

    val results = receiveN(2, 2 seconds).toSet
    val mixedResponses = results collect {
      case msg: DockerInfoSuccessResponse => msg
      // Scoop out the request here since we can't match the exception on the whole message.
      case msg: WorkflowDockerLookupFailure if msg.reason.getMessage == "Failed to get docker hash for ubuntu:older Lookup failed" => msg.request
    }

    Set(LatestSuccessResponse, OlderRequest) should equal(mixedResponses)

    // Try again, I have a good feeling about this.
    lookupActor ! OlderRequest
    dockerHashingActor.expectMsg(OlderRequest)
    dockerHashingActor.reply(OlderSuccessResponse)
    expectMsg(OlderSuccessResponse)
  }

  it should "reuse previously looked up hashes following a restart" in {
    val db = dbWithQuery {
      Future.successful(
        Seq(LatestStoreEntry(workflowId), OlderStoreEntry(workflowId)))
    }

    val lookupActor = TestActorRef(WorkflowDockerLookupActor.props(workflowId, dockerHashingActor.ref, isRestart = true, db))

    lookupActor ! LatestRequest
    lookupActor ! OlderRequest

    dockerHashingActor.expectNoMessage()

    val results = receiveN(2, 2 seconds).toSet
    val successes = results collect { case result: DockerInfoSuccessResponse => result }

    successes should equal(Set(LatestSuccessResponse, OlderSuccessResponse))
  }

  it should "not try to look up hashes if not restarting" in {
    val db = dbWithWrite(Future.successful(()))
    val lookupActor = TestActorRef(WorkflowDockerLookupActor.props(workflowId, dockerHashingActor.ref, isRestart = false, db))

    lookupActor ! LatestRequest
    lookupActor ! OlderRequest

    dockerHashingActor.expectMsg(LatestRequest)
    dockerHashingActor.expectMsg(OlderRequest)
    dockerHashingActor.reply(LatestSuccessResponse)
    dockerHashingActor.reply(OlderSuccessResponse)

    val results = receiveN(2, 2 seconds).toSet
    val successes = results collect { case result: DockerInfoSuccessResponse => result }

    successes should equal(Set(LatestSuccessResponse, OlderSuccessResponse))
  }

  it should "handle hash write errors appropriately" in {
    val db = dbWithWrite {
      numWrites = numWrites + 1
      if (numWrites == 1) Future.failed(new RuntimeException("Fake exception from a test.")) else Future.successful(())
    }

    val lookupActor = TestActorRef(WorkflowDockerLookupActor.props(workflowId, dockerHashingActor.ref, isRestart = false, db))
    lookupActor ! LatestRequest

    // The WorkflowDockerLookupActor should not have the hash for this tag yet and will need to query the dockerHashingActor.
    dockerHashingActor.expectMsg(LatestRequest)
    dockerHashingActor.reply(LatestSuccessResponse)
    // The WorkflowDockerLookupActor is going to fail when it tries to write to that broken DB.
    expectMsgClass(classOf[WorkflowDockerLookupFailure])
    numWrites should equal(1)

    lookupActor ! LatestRequest
    // The WorkflowDockerLookupActor will query the dockerHashingActor again.
    dockerHashingActor.expectMsg(LatestRequest)
    dockerHashingActor.reply(LatestSuccessResponse)
    // The WorkflowDockerLookupActor should forward the success message to this actor.
    expectMsg(LatestSuccessResponse)
    numWrites should equal(2)
  }

  it should "emit a terminal failure message if failing to read hashes on restart" in {
    val db = dbWithQuery {
      numReads = numReads + 1
      Future.failed(new Exception("Don't worry this is just a dummy failure in a test") with NoStackTrace)
    }

    val lookupActor = TestActorRef(WorkflowDockerLookupActor.props(workflowId, dockerHashingActor.ref, isRestart = true, db))
    lookupActor ! LatestRequest

    dockerHashingActor.expectNoMessage()
    expectMsgClass(classOf[WorkflowDockerTerminalFailure])
    numReads should equal(1)
  }

  it should "emit a terminal failure message if unable to parse hashes read from the database on restart" in {
    val db = dbWithQuery {
      numReads = numReads + 1
      Future.successful(Seq(
        DockerHashStoreEntry(workflowId.toString, Latest, "md5:AAAAA", None),
        // missing the "algorithm:" preceding the hash value so this should fail parsing.
        DockerHashStoreEntry(workflowId.toString, Older, "BBBBB", None)
      ))
    }

    val lookupActor = TestActorRef(WorkflowDockerLookupActor.props(workflowId, dockerHashingActor.ref, isRestart = true, db))
    lookupActor ! LatestRequest

    dockerHashingActor.expectNoMessage()
    expectMsgClass(classOf[WorkflowDockerTerminalFailure])
    numReads should equal(1)
  }

  def dbWithWrite(writeFn: => Future[Unit]): EngineSlickDatabase = {
    databaseInterface(write = _ => writeFn)
  }

  def dbWithQuery(queryFn: => Future[Seq[DockerHashStoreEntry]]): EngineSlickDatabase = {
    databaseInterface(query = _ => queryFn)
  }

  def databaseInterface(query: String => Future[Seq[DockerHashStoreEntry]] = abjectFailure,
                        write: DockerHashStoreEntry => Future[Unit] = abjectFailure): EngineSlickDatabase = {
    new EngineSlickDatabase(DatabaseConfig) {
      override def queryDockerHashStoreEntries(workflowExecutionUuid: String)(implicit ec: ExecutionContext): Future[Seq[DockerHashStoreEntry]] = query(workflowExecutionUuid)

      override def addDockerHashStoreEntry(dockerHashStoreEntry: DockerHashStoreEntry)(implicit ec: ExecutionContext): Future[Unit] = write(dockerHashStoreEntry)
    }.initialized(EngineServicesStore.EngineLiquibaseSettings)
  }
}


object WorkflowDockerLookupActorSpec {
  val Latest = "ubuntu:latest"
  val Older = "ubuntu:older"

  val LatestImageId = DockerImageIdentifier.fromString(Latest).get.asInstanceOf[DockerImageIdentifierWithoutHash]
  val OlderImageId = DockerImageIdentifier.fromString(Older).get.asInstanceOf[DockerImageIdentifierWithoutHash]

  val LatestRequest = DockerInfoRequest(LatestImageId)
  val OlderRequest = DockerInfoRequest(OlderImageId)

  def LatestStoreEntry(workflowId: WorkflowId): DockerHashStoreEntry = DockerHashStoreEntry(workflowId.toString, Latest, "md5:AAAAAAAA", None)
  def OlderStoreEntry(workflowId: WorkflowId): DockerHashStoreEntry = DockerHashStoreEntry(workflowId.toString, Older, "md5:BBBBBBBB", None)

  val LatestSuccessResponse = DockerInfoSuccessResponse(DockerInformation(DockerHashResult("md5", "AAAAAAAA"), None), LatestRequest)
  val OlderSuccessResponse = DockerInfoSuccessResponse(DockerInformation(DockerHashResult("md5", "BBBBBBBB"), None), OlderRequest)

  val DatabaseConfig = ConfigFactory.load.getConfig("database")

  def abjectFailure[A, B]: A => Future[B] = _ => Future.failed(new RuntimeException("Should not be called!"))

  class TestWorkflowDockerLookupActor(workflowId: WorkflowId, dockerHashingActor: ActorRef, startState: StartableState, backoff: Backoff)
    extends WorkflowDockerLookupActor(
      workflowId,
      dockerHashingActor,
      startState.restarted,
      EngineServicesStore.engineDatabaseInterface) {
    override protected def initialBackoff = backoff
  }
}
