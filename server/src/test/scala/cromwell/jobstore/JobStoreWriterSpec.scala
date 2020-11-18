package cromwell.jobstore

import akka.testkit.{TestFSMRef, TestProbe}
import common.assertion.CromwellTimeoutSpec
import common.collections.WeightedQueue
import cromwell.CromwellTestKitWordSpec
import cromwell.core.actor.BatchActor.{BatchActorState, CommandAndReplyTo, Processing}
import cromwell.core.{CallOutputs, WorkflowId}
import cromwell.engine.workflow.{CoordinatedWorkflowStoreActorBuilder, SqlWorkflowStoreBuilder}
import cromwell.jobstore.JobStore.{JobCompletion, WorkflowCompletion}
import cromwell.jobstore.JobStoreActor.{JobStoreWriteSuccess, JobStoreWriterCommand, RegisterJobCompleted, RegisterWorkflowCompleted}
import org.scalatest.BeforeAndAfter
import org.scalatest.matchers.should.Matchers
import wom.graph.GraphNodePort.OutputPort

import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, Future, Promise}
import scala.language.postfixOps

class JobStoreWriterSpec extends CromwellTestKitWordSpec with SqlWorkflowStoreBuilder with CoordinatedWorkflowStoreActorBuilder with Matchers with BeforeAndAfter with CromwellTimeoutSpec {

  var database: WriteCountingJobStore = _
  var workflowId: WorkflowId = _
  val successResult: JobResult = JobResultSuccess(Some(0), CallOutputs.empty)
  private val flushFrequency = 0.5.seconds

  before {
    database = WriteCountingJobStore.makeNew
    workflowId = WorkflowId.randomId()
  }

  private def sendRegisterCompletion(jobStoreWriter: TestFSMRef[BatchActorState, WeightedQueue[CommandAndReplyTo[JobStoreWriterCommand], Int], JobStoreWriterActor])(attempt: Int): Unit = {
    jobStoreWriter ! RegisterJobCompleted(jobKey(attempt), successResult)
  }

  private def jobKey(attempt: Int): JobStoreKey = JobStoreKey(workflowId, s"call.fqn", None, attempt)

  private def assertWriteSuccess(key: JobStoreKey, result: JobResult): Unit = {
    key.workflowId shouldBe workflowId
    key.callFqn shouldBe "call.fqn"
    key.index shouldBe None
    result shouldBe successResult
    ()
  }

  //noinspection SameParameterValue
  private def assertDb(totalWritesCalled: Int, jobCompletionsRecorded: Int, workflowCompletionsRecorded: Int): Unit = {
    database.totalWritesCalled shouldBe totalWritesCalled
    database.jobCompletionsRecorded shouldBe jobCompletionsRecorded
    database.workflowCompletionsRecorded shouldBe workflowCompletionsRecorded
    ()
  }

  //noinspection SameParameterValue
  private def assertReceived(expectedJobStoreWriteAcks: Int): Unit = {
    val received = receiveN(expectedJobStoreWriteAcks, 10 seconds)
    received foreach {
      case JobStoreWriteSuccess(RegisterJobCompleted(key: JobStoreKey, result)) => assertWriteSuccess(key, result)
      case JobStoreWriteSuccess(RegisterWorkflowCompleted(id)) => id shouldBe workflowId
      case message => fail(s"Unexpected response message: $message")
    }
    ()
  }

  "JobStoreWriter" should {
    "be able to collapse writes together if they arrive while a database access is ongoing" in {
      runWithDatabase(databaseConfig) { workflowStore =>
        val jobStoreWriter =
          TestFSMRef(
            factory = new JobStoreWriterActor(
              jsd = database,
              batchSize = 5,
              flushRate = flushFrequency,
              serviceRegistryActor = TestProbe("serviceRegistryActor-collapse").ref,
              threshold = 1000,
              workflowStoreAccess = access("coordinatedAccessActor-collapse")(workflowStore),
            ),
            name = "jobStoreWriter-collapse",
          )

        // Send a job completion. The database will hang.
        sendRegisterCompletion(jobStoreWriter)(1)
        val writer = jobStoreWriter.underlyingActor
        eventually {
          writer.stateName shouldBe Processing
        }

        // Send some more completions.  These should pile up in the state data.
        List(2, 3) foreach sendRegisterCompletion(jobStoreWriter)

        writer.stateData.weight should equal(2)
        database.continue()

        assertReceived(expectedJobStoreWriteAcks = 3)

        assertDb(
          totalWritesCalled = 2,
          jobCompletionsRecorded = 3,
          workflowCompletionsRecorded = 0
        )
      }
    }

    "be able to skip job-completion writes if the workflow completes, but still respond appropriately" in {
      runWithDatabase(databaseConfig) { workflowStore =>
        val jobStoreWriter =
          TestFSMRef(
            new JobStoreWriterActor(
              jsd = database,
              batchSize = 5,
              flushRate = flushFrequency,
              serviceRegistryActor = TestProbe("serviceRegistryActor-skip").ref,
              threshold = 1000,
              workflowStoreAccess = access("coordinatedAccessActor-skip")(workflowStore),
            ),
            "jobStoreWriter-skip",
          )

        // Send a job completion. The database will hang.
        sendRegisterCompletion(jobStoreWriter)(1)
        val writer = jobStoreWriter.underlyingActor
        eventually {
          writer.stateName shouldBe Processing
        }

        // Send some more completions.  These should pile up in the state data.
        List(2, 3) foreach sendRegisterCompletion(jobStoreWriter)
        jobStoreWriter ! RegisterWorkflowCompleted(workflowId)

        writer.stateData.weight should equal(3)
        // The testing DB intentionally blocks after the first write until `continue` is called.
        database.continue()

        assertReceived(expectedJobStoreWriteAcks = 3)

        assertDb(
          totalWritesCalled = 2,
          jobCompletionsRecorded = 1,
          workflowCompletionsRecorded = 1
        )
      }
    }
  }
}

class WriteCountingJobStore(var totalWritesCalled: Int, var jobCompletionsRecorded: Int, var workflowCompletionsRecorded: Int) extends JobStore {

  // A Promise so that the calling tests can hang the writer on the db write.  Once the promise is completed the writer is
  // released and all further messages will be written immediately.
  private val writePromise = Promise[Unit]()

  def continue(): Unit = {
    writePromise.trySuccess(())
    ()
  }

  override def writeToDatabase(workflowCompletions: Seq[WorkflowCompletion], jobCompletions: Seq[JobCompletion], batchSize: Int)
                              (implicit ec: ExecutionContext): Future[Unit] = {

    totalWritesCalled += 1
    jobCompletionsRecorded += jobCompletions.size
    workflowCompletionsRecorded += workflowCompletions.size
    writePromise.future
  }

  override def readJobResult(jobStoreKey: JobStoreKey, taskOutputs: Seq[OutputPort])(implicit ec: ExecutionContext): Future[Option[JobResult]] = throw new UnsupportedOperationException()
}

object WriteCountingJobStore {
  def makeNew = new WriteCountingJobStore(0, 0, 0)
}
