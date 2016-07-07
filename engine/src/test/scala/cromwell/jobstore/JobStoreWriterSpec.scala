package cromwell.jobstore

import cromwell.CromwellTestkitSpec
import cromwell.core.WorkflowId
import org.scalatest.Matchers
import scala.concurrent.duration._
import language.postfixOps

import scala.concurrent.Future

class JobStoreWriterSpec extends CromwellTestkitSpec with Matchers {

  "JobStoreWriter" should {
    "be able to collapse writes together if they arrive while a database access is ongoing" in {
      val database = WriteCountingJobStoreDatabase.makeNew
      val jsw = system.actorOf(JobStoreWriter.props(database))
      val wfid = WorkflowId.randomId()
      def jobKey(number: Int): JobStoreKey = JobStoreKey(wfid, s"call.fqn_$number", None)
      val jobResult: JobResult = JobResultSuccess(Some(0), Map.empty)

      jsw ! RegisterJobCompleted(jobKey(0), jobResult)
      jsw ! RegisterJobCompleted(jobKey(1), jobResult)
      jsw ! RegisterJobCompleted(jobKey(2), jobResult)
      val received = receiveN(3, 10 seconds)
      received foreach {
        case JobStoreWriteSuccess(RegisterJobCompleted(JobStoreKey(jsk_wfid, jsk_callfqn, None), result)) =>
          jsk_wfid shouldBe wfid
          jsk_callfqn should startWith("call.fqn")
          result shouldBe jobResult
        case message => fail(s"Unexpected response message: $message")
      }

      database.totalWritesCalled shouldBe 2
      database.jobCompletionsRecorded shouldBe 3
      database.workflowCompletionsRecorded shouldBe 0
    }

    "be able to skip job-completion writes if the workflow completes, but still respond appropriately" in {
      val database = WriteCountingJobStoreDatabase.makeNew
      val jsw = system.actorOf(JobStoreWriter.props(database))
      val wfid = WorkflowId.randomId()
      def jobKey(number: Int): JobStoreKey = JobStoreKey(wfid, s"call.fqn_$number", None)
      val jobResult: JobResult = JobResultSuccess(Some(0), Map.empty)

      jsw ! RegisterJobCompleted(jobKey(0), jobResult)
      jsw ! RegisterJobCompleted(jobKey(1), jobResult)
      jsw ! RegisterWorkflowCompleted(wfid)
      val received = receiveN(3, 10 seconds)
      received foreach {
        case JobStoreWriteSuccess(RegisterJobCompleted(JobStoreKey(jsk_wfid, jsk_callfqn, None), result)) =>
          jsk_wfid shouldBe wfid
          jsk_callfqn should startWith("call.fqn_")
          result shouldBe jobResult
        case JobStoreWriteSuccess(RegisterWorkflowCompleted(rwc_wfid)) =>
          rwc_wfid shouldBe wfid
        case message => fail(s"Unexpected response message: $message")
      }

      database.totalWritesCalled shouldBe 2
      database.jobCompletionsRecorded shouldBe 1
      database.workflowCompletionsRecorded shouldBe 1
    }
  }
}

class WriteCountingJobStoreDatabase(var totalWritesCalled: Int, var jobCompletionsRecorded: Int, var workflowCompletionsRecorded: Int) extends JobStoreDatabase {

  override def writeToDatabase(jobCompletions: Map[JobStoreKey, JobResult], workflowCompletions: List[WorkflowId]): Future[Unit] = {
    totalWritesCalled += 1
    jobCompletionsRecorded += jobCompletions.size
    workflowCompletionsRecorded += workflowCompletions.size
    Future.successful(())
  }
}

object WriteCountingJobStoreDatabase {
  def makeNew = new WriteCountingJobStoreDatabase(0, 0, 0)
}