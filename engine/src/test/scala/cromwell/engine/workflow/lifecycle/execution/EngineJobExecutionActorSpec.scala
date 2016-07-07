package cromwell.engine.workflow.lifecycle.execution

import akka.actor.{Actor, Props}
import akka.testkit.{TestFSMRef, TestProbe}
import cromwell.CromwellTestkitSpec
import cromwell.backend.BackendJobExecutionActor._
import cromwell.backend._
import cromwell.core.{ExecutionStore, JobOutputs, OutputStore, WorkflowId}
import cromwell.engine.workflow.WorkflowDescriptorBuilder
import cromwell.engine.workflow.lifecycle.execution.EngineJobExecutionActor._
import cromwell.engine.workflow.lifecycle.execution.JobPreparationActor.BackendJobPreparationFailed
import cromwell.jobstore.{JobComplete, JobNotComplete, JobResultFailure, JobResultSuccess}
import cromwell.util.SampleWdl
import org.scalatest.{BeforeAndAfterAll, Matchers}
import org.specs2.mock.Mockito
import wdl4s.expression.{NoFunctions, WdlStandardLibraryFunctions}
import wdl4s.{Call, Task}

class EngineJobExecutionActorSpec extends CromwellTestkitSpec with Matchers with WorkflowDescriptorBuilder with Mockito with BeforeAndAfterAll {

  override implicit val actorSystem = system

  private val workflowId: WorkflowId = WorkflowId.randomId()
  val descriptor = createMaterializedEngineWorkflowDescriptor(workflowId, SampleWdl.HelloWorld.asWorkflowSources())
  val backendProbe = TestProbe()
  val mockBackendProps = Props(new Actor {
    def receive = {
      case x => backendProbe.ref forward x
    }
  })
  val ejeaParent = TestProbe()
  val factory = new BackendLifecycleActorFactory {
    override def workflowInitializationActorProps(workflowDescriptor: BackendWorkflowDescriptor, calls: Seq[Call]): Option[Props] = None
    override def jobExecutionActorProps(jobDescriptor: BackendJobDescriptor, initializationData: Option[BackendInitializationData]): Props = mockBackendProps
    override def expressionLanguageFunctions(workflowDescriptor: BackendWorkflowDescriptor, jobKey: BackendJobDescriptorKey, initializationData: Option[BackendInitializationData]): WdlStandardLibraryFunctions = {
      NoFunctions
    }
  }

  def buildEJEA(restarting: Boolean) = new TestFSMRef[EngineJobExecutionActorState, Unit, EngineJobExecutionActor](system, EngineJobExecutionActor.props(
    WorkflowExecutionActorData(descriptor, ExecutionStore(Map.empty), Map.empty, OutputStore(Map.empty)),
    factory,
    None,
    restarting = restarting
  ), ejeaParent.ref, s"EngineJobExecutionActorSpec-$workflowId")

  "EngineJobExecutionActorSpec" should {
    "send a Job SucceededResponse if the job is already complete and successful" in {
      val ejea = buildEJEA(restarting = true)
      ejea.setState(CheckingJobStatus)
      val returnCode: Option[Int] = Option(0)
      val jobOutputs: JobOutputs = Map.empty

      ejea ! JobComplete(mock[BackendJobDescriptorKey], JobResultSuccess(returnCode, jobOutputs))

      ejeaParent.expectMsgPF(awaitTimeout) {
        case response: SucceededResponse =>
          response.returnCode shouldBe returnCode
          response.jobOutputs shouldBe jobOutputs
      }

      ejea.stop()
    }

    "send a Job FailedNonRetryableResponse if the job is already complete and failed" in {
      val ejea = buildEJEA(restarting = true)
      ejea.setState(CheckingJobStatus)
      val returnCode: Option[Int] = Option(1)
      val reason: Throwable = new Exception("something horrible happened...")

      ejea ! JobComplete(mock[BackendJobDescriptorKey], JobResultFailure(returnCode, reason))

      ejeaParent.expectMsgPF(awaitTimeout) {
        case response: FailedNonRetryableResponse =>
          response.returnCode shouldBe returnCode
          response.throwable shouldBe reason
      }

      ejea.stop()
    }

    "send a RecoverJobCommand to the backend if the job is not complete and EJEA is in restarting mode" in {
      val ejea = buildEJEA(restarting = true)
      ejea.setState(CheckingJobStatus)
      val task = mock[Task]
      task.declarations returns Seq.empty

      val jobKey = BackendJobDescriptorKey(Call(None, "wf.call", task, Set.empty, Map.empty, None), None, 1)
      ejea ! JobNotComplete(jobKey)

      backendProbe.expectMsg(awaitTimeout, RecoverJobCommand)

      ejea.stop()
    }

    // This test is bit redundant with the one above, hoever this one tests the behaviour when getting a Restart message
    // whereas the test above tests the behaviour when getting a JobNotComplete response from JobReader
    "send an RecoverJobCommand to the backend when getting an Execute message if the Job is not complete and EJEA is in restarting mode" in {
      val ejea = buildEJEA(restarting = true)
      ejea.setState(Pending)
      val task = mock[Task]
      task.declarations returns Seq.empty

      val jobKey = BackendJobDescriptorKey(Call(None, "wf.call", task, Set.empty, Map.empty, None), None, 1)
      ejea ! EngineJobExecutionActor.Execute(jobKey)

      backendProbe.expectMsg(awaitTimeout, RecoverJobCommand)

      ejea.stop()
    }

    "send an ExecuteJobCommand to the backend when getting an Execute message if the Job is not complete and EJEA is NOT in restarting mode" in {
      val ejea = buildEJEA(restarting = false)
      ejea.setState(Pending)
      val task = mock[Task]
      task.declarations returns Seq.empty

      val jobKey = BackendJobDescriptorKey(Call(None, "wf.call", task, Set.empty, Map.empty, None), None, 1)
      ejea ! EngineJobExecutionActor.Execute(jobKey)

      backendProbe.expectMsg(awaitTimeout, ExecuteJobCommand)

      ejea.stop()
    }

    "forward JobPreparationFailed to its parent and die" in {
      val deathWatch = TestProbe()
      val jobKey = mock[BackendJobDescriptorKey]

      val ejea = buildEJEA(restarting = true)
      ejea.setState(PreparingJob)

      val failed = BackendJobPreparationFailed(jobKey, new scala.Exception("Failed"))
      deathWatch watch ejea

      ejea ! failed
      ejeaParent.expectMsg(awaitTimeout, failed)
      ejeaParent.lastSender shouldBe self

      deathWatch.expectTerminated(ejea)
    }

    "forward responds from backend to its parent" in {
      val deathWatch = TestProbe()

      val jobKey = mock[BackendJobDescriptorKey]

      def ejeaInRunningState() = {
        val ejea = buildEJEA(restarting = true)
        ejea.setState(RunningJob)
        ejea
      }

      val successResponse = SucceededResponse(jobKey, Option(0), Map.empty)
      val failureRetryableResponse = FailedRetryableResponse(jobKey, new Exception("Failed"), Option(0))
      val failureNonRetryableResponse = FailedNonRetryableResponse(jobKey, new Exception("Failed"), Option(0))
      val abortedResponse = FailedNonRetryableResponse(jobKey, new Exception("Failed"), Option(0))

      val ejea1 = ejeaInRunningState()
      deathWatch watch ejea1
      ejea1 ! successResponse
      ejeaParent.expectMsg(awaitTimeout, successResponse)
      ejeaParent.lastSender shouldBe self
      deathWatch.expectTerminated(ejea1)

      val ejea2 = ejeaInRunningState()
      deathWatch watch ejea2
      ejea2 ! failureRetryableResponse
      ejeaParent.expectMsg(awaitTimeout, failureRetryableResponse)
      ejeaParent.lastSender shouldBe self
      deathWatch.expectTerminated(ejea2)

      val ejea3 = ejeaInRunningState()
      deathWatch watch ejea3
      ejea3 ! failureNonRetryableResponse
      ejeaParent.expectMsg(awaitTimeout, failureNonRetryableResponse)
      ejeaParent.lastSender shouldBe self
      deathWatch.expectTerminated(ejea3)

      val ejea4 = ejeaInRunningState()
      deathWatch watch ejea4
      ejea4 ! abortedResponse
      ejeaParent.expectMsg(awaitTimeout, abortedResponse)
      ejeaParent.lastSender shouldBe self
      deathWatch.expectTerminated(ejea4)
    }
  }

  override def afterAll(): Unit = {
    system.shutdown()
  }

}
