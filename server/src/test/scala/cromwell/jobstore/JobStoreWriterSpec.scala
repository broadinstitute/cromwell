package cromwell.jobstore

import akka.testkit.{TestFSMRef, TestProbe}
import common.collections.WeightedQueue
import cromwell.CromwellTestKitWordSpec
import cromwell.core.actor.BatchActor.{BatchActorState, CommandAndReplyTo, Processing}
import cromwell.core.{CallOutputs, WorkflowId}
import cromwell.jobstore.JobStore.{JobCompletion, WorkflowCompletion}
import cromwell.jobstore.JobStoreActor.{JobStoreWriteSuccess, JobStoreWriterCommand, RegisterJobCompleted, RegisterWorkflowCompleted}
import org.scalatest.{BeforeAndAfter, Matchers}
import wom.graph.GraphNodePort.OutputPort

import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, Future, Promise}
import scala.language.postfixOps

class JobStoreWriterSpec extends CromwellTestKitWordSpec with Matchers with BeforeAndAfter {
  
  var database: WriteCountingJobStore = _
  var jobStoreWriter: TestFSMRef[BatchActorState, WeightedQueue[CommandAndReplyTo[JobStoreWriterCommand], Int], JobStoreWriterActor] = _
  var workflowId: WorkflowId = _
  val successResult: JobResult = JobResultSuccess(Some(0), CallOutputs.empty)
  val flushFrequency = 0.5 second
  val sleepMillis = flushFrequency.toMillis * 5

  before {
    database = WriteCountingJobStore.makeNew
    jobStoreWriter = TestFSMRef(new JobStoreWriterActor(database, 5, flushFrequency, TestProbe().ref, 1000))
    workflowId = WorkflowId.randomId()
  }

  private def sendRegisterCompletion(attempt: Int): Unit = {
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

  private def assertDb(totalWritesCalled: Int, jobCompletionsRecorded: Int, workflowCompletionsRecorded: Int): Unit = {
    database.totalWritesCalled shouldBe totalWritesCalled
    database.jobCompletionsRecorded shouldBe jobCompletionsRecorded
    database.workflowCompletionsRecorded shouldBe workflowCompletionsRecorded
    ()
  }

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

      // Send a job completion. The database will hang.
      sendRegisterCompletion(1)
      val writer = jobStoreWriter.underlyingActor
      eventually {
        writer.stateName shouldBe Processing
      }

      // Send some more completions.  These should pile up in the state data.
      List(2, 3) foreach sendRegisterCompletion

      writer.stateData.weight should equal(2)
      database.continue()

      assertReceived(expectedJobStoreWriteAcks = 3)

      assertDb(
        totalWritesCalled = 2,
        jobCompletionsRecorded = 3,
        workflowCompletionsRecorded = 0
      )
    }

    "be able to skip job-completion writes if the workflow completes, but still respond appropriately" in {

      // Send a job completion. The database will hang.
      sendRegisterCompletion(1)
      val writer = jobStoreWriter.underlyingActor
      eventually {
        writer.stateName shouldBe Processing
      }

      // Send some more completions.  These should pile up in the state data.
      List(2, 3) foreach sendRegisterCompletion
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

class WriteCountingJobStore(var totalWritesCalled: Int, var jobCompletionsRecorded: Int, var workflowCompletionsRecorded: Int) extends JobStore {

  // A Promise so that the calling tests can hang the writer on the db write.  Once the promise is completed the writer is
  // released and all further messages will be written immediately.
  val writePromise = Promise[Unit]()

  def continue() = writePromise.trySuccess(())

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
