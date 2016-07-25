package cromwell.jobstore

import akka.actor.Props
import com.typesafe.config.ConfigFactory
import cromwell.CromwellTestkitSpec
import cromwell.backend.BackendJobDescriptorKey
import cromwell.core.{JobOutput, WorkflowId}
import cromwell.jobstore.JobStoreActor._
import cromwell.jobstore.JobStoreServiceSpec._
import org.scalatest.Matchers
import org.specs2.mock.Mockito
import wdl4s.Call
import wdl4s.values.WdlString

import scala.concurrent.duration._
import scala.language.postfixOps

object JobStoreServiceSpec {
  val MaxWait = 5 seconds
}

class JobStoreServiceSpec extends CromwellTestkitSpec with Matchers with Mockito {

  "JobStoreService" should {
    "work" in {
      val config = ConfigFactory.parseString("{}")
      val jobStoreService = system.actorOf(JobStoreActor.props(WriteCountingJobStoreDatabase.makeNew))

      val workflowId = WorkflowId.randomId()
      val successCall = mock[Call]
      successCall.fullyQualifiedName returns "foo.bar"
      val successKey = BackendJobDescriptorKey(successCall, None, 1).toJobStoreKey(workflowId)

      jobStoreService ! QueryJobCompletion(successKey)
      expectMsgType[JobNotComplete.type](MaxWait)

      val outputs = Map("baz" -> JobOutput(WdlString("qux")))

      jobStoreService ! RegisterJobCompleted(successKey, JobResultSuccess(Option(0), outputs))
      expectMsgType[JobStoreWriteSuccess](MaxWait)

      jobStoreService ! QueryJobCompletion(successKey)
      expectMsgPF(MaxWait) {
        case JobComplete(JobResultSuccess(Some(0), os)) if os == outputs =>
      }

      val failureCall = mock[Call]
      failureCall.fullyQualifiedName returns "baz.qux"
      val failureKey = BackendJobDescriptorKey(failureCall, None, 1).toJobStoreKey(workflowId)

      jobStoreService ! QueryJobCompletion(failureKey)
      expectMsgType[JobNotComplete.type](MaxWait)

      jobStoreService ! RegisterJobCompleted(failureKey, JobResultFailure(Option(11), new IllegalArgumentException("Insufficient funds")))
      expectMsgType[JobStoreWriteSuccess](MaxWait)

      jobStoreService ! QueryJobCompletion(failureKey)
      expectMsgPF(MaxWait) {
        case JobComplete(JobResultFailure(Some(11), _)) =>
      }

      jobStoreService ! RegisterWorkflowCompleted(workflowId)
      expectMsgType[JobStoreWriteSuccess](MaxWait)
    }
  }
}
