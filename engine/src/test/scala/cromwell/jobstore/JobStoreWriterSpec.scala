package cromwell.jobstore

import akka.testkit.TestFSMRef
import cromwell.CromwellTestKitSpec
import cromwell.core.WorkflowId
import cromwell.jobstore.JobStoreActor.{JobStoreWriteSuccess, RegisterJobCompleted, RegisterWorkflowCompleted}
import org.scalatest.{BeforeAndAfter, Matchers}
import wdl4s.TaskOutput

import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, Future, Promise}
import scala.language.postfixOps

class JobStoreWriterSpec extends CromwellTestKitSpec with Matchers with BeforeAndAfter {
  
  var database: WriteCountingJobStore = _
  var jobStoreWriter: TestFSMRef[JobStoreWriterState, JobStoreWriterData, JobStoreWriterActor] = _
  var workflowId: WorkflowId = _
  val successResult: JobResult = JobResultSuccess(Some(0), Map.empty)

  before {
    database = WriteCountingJobStore.makeNew
    jobStoreWriter = TestFSMRef(new JobStoreWriterActor(database))
    workflowId = WorkflowId.randomId()
  }

  private def sendRegisterCompletions(attempts: Int): Unit = {
    0 until attempts foreach { a => jobStoreWriter ! RegisterJobCompleted(jobKey(attempt = a), successResult) }
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
    jobStoreWriter.underlyingActor.stateName shouldBe Pending
    ()
  }

  "JobStoreWriter" should {
    "be able to collapse writes together if they arrive while a database access is ongoing" in {

      sendRegisterCompletions(attempts = 3)

      val writer = jobStoreWriter.underlyingActor
      writer.stateName shouldBe WritingToDatabase
      writer.stateData.currentOperation should have size 1
      writer.stateData.nextOperation should have size 2

      // The testing DB intentionally blocks after the first write until `continue` is called.
      database.continue()

      assertReceived(expectedJobStoreWriteAcks = 3)
      writer.stateData.currentOperation shouldBe empty
      writer.stateData.nextOperation shouldBe empty

      assertDb(
        totalWritesCalled = 2,
        jobCompletionsRecorded = 3,
        workflowCompletionsRecorded = 0
      )
    }

    "be able to skip job-completion writes if the workflow completes, but still respond appropriately" in {

      sendRegisterCompletions(attempts = 2)
      jobStoreWriter ! RegisterWorkflowCompleted(workflowId)

      val writer = jobStoreWriter.underlyingActor
      writer.stateName shouldBe WritingToDatabase
      writer.stateData.currentOperation should have size 1
      writer.stateData.nextOperation should have size 2

      // The testing DB intentionally blocks after the first write until `continue` is called.
      database.continue()

      assertReceived(expectedJobStoreWriteAcks = 3)
      writer.stateData.currentOperation shouldBe empty
      writer.stateData.nextOperation shouldBe empty

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

  override def writeToDatabase(jobCompletions: Map[JobStoreKey, JobResult], workflowCompletions: List[WorkflowId])
                              (implicit ec: ExecutionContext): Future[Unit] = {
    totalWritesCalled += 1
    jobCompletionsRecorded += jobCompletions.size
    workflowCompletionsRecorded += workflowCompletions.size
    writePromise.future
  }

  override def readJobResult(jobStoreKey: JobStoreKey, taskOutputs: Seq[TaskOutput])(implicit ec: ExecutionContext): Future[Option[JobResult]] = throw new NotImplementedError()
}

object WriteCountingJobStore {
  def makeNew = new WriteCountingJobStore(0, 0, 0)
}
