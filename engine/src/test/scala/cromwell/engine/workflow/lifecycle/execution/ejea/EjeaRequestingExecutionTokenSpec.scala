package cromwell.engine.workflow.lifecycle.execution.ejea

import cromwell.engine.workflow.lifecycle.execution.EngineJobExecutionActor._
import cromwell.engine.workflow.lifecycle.execution.JobStarting
import cromwell.engine.workflow.lifecycle.execution.WorkflowExecutionActor.RequestOutputStore
import cromwell.engine.workflow.tokens.JobExecutionTokenDispenserActor.{JobExecutionTokenDenied, JobExecutionTokenDispensed}
import cromwell.jobstore.JobStoreActor.QueryJobCompletion
import org.scalatest.concurrent.Eventually

class EjeaRequestingExecutionTokenSpec extends EngineJobExecutionActorSpec with CanValidateJobStoreKey with Eventually {

  override implicit val stateUnderTest: EngineJobExecutionActorState = RequestingExecutionToken

  "An EJEA in the RequestingExecutionToken state" should {

    List(true, false) foreach { restarting =>
      s"do nothing when denied a token (with restarting=$restarting)" in {
        ejea = helper.buildEJEA(restarting = restarting)
        ejea ! JobExecutionTokenDenied(1) // 1 is arbitrary. Doesn't matter what position in the queue we are.

        helper.jobTokenDispenserProbe.expectNoMsg(max = awaitAlmostNothing)
        helper.jobPreparationProbe.msgAvailable should be(false)
        helper.jobStoreProbe.msgAvailable should be(false)

        ejea.stateName should be(RequestingExecutionToken)
      }
    }

    CallCachingModes foreach { mode =>
      s"check against the Job Store if restarting is true ($mode)" in {
        ejea = helper.buildEJEA(restarting = true)
        ejea ! JobExecutionTokenDispensed(helper.executionToken)

        helper.jobStoreProbe.expectMsgPF(max = awaitTimeout, hint = "Awaiting job store lookup") {
          case QueryJobCompletion(jobKey, taskOutputs) =>
            validateJobStoreKey(jobKey)
            taskOutputs should be(helper.task.outputs)
        }
        helper.bjeaProbe.expectNoMsg(awaitAlmostNothing)
        helper.jobHashingInitializations shouldBe NothingYet
        ejea.stateName should be(CheckingJobStore)
      }


      s"bypass the Job Store and request output store to start preparing the job for running or call caching ($mode)" in {
        ejea = helper.buildEJEA(restarting = false)
        ejea ! JobExecutionTokenDispensed(helper.executionToken)

        helper.replyToProbe.expectMsg(max = awaitTimeout, hint = "Awaiting JobStarting message", JobStarting(helper.jobDescriptorKey))
        helper.replyToProbe.expectMsg(max = awaitTimeout, hint = "Awaiting RequestOutputStore message", RequestOutputStore)
        ejea.stateName should be(WaitingForOutputStore)
      }
    }
  }
}
