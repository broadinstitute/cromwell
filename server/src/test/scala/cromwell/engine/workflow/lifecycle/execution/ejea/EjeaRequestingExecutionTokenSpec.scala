package cromwell.engine.workflow.lifecycle.execution.ejea

import cromwell.engine.workflow.lifecycle.execution.JobStarting
import cromwell.engine.workflow.lifecycle.execution.WorkflowExecutionActor.RequestValueStore
import cromwell.engine.workflow.lifecycle.execution.job.EngineJobExecutionActor._
import cromwell.engine.workflow.tokens.JobTokenDispenserActor.JobTokenDispensed
import cromwell.jobstore.JobStoreActor.QueryJobCompletion
import org.scalatest.concurrent.Eventually

class EjeaRequestingExecutionTokenSpec extends EngineJobExecutionActorSpec with CanValidateJobStoreKey with Eventually {

  override implicit val stateUnderTest: EngineJobExecutionActorState = RequestingExecutionToken

  "An EJEA in the RequestingExecutionToken state" should {

    List(true, false) foreach { restarting =>
      s"do nothing when denied a token (with restarting=$restarting)" in {
        ejea = helper.buildEJEA(restarting = restarting)

        helper.jobExecutionTokenDispenserProbe.expectNoMessage(max = awaitAlmostNothing)
        helper.jobPreparationProbe.msgAvailable should be(false)
        helper.jobStoreProbe.msgAvailable should be(false)

        ejea.stateName should be(RequestingExecutionToken)
      }
    }

    CallCachingModes foreach { mode =>
      s"check against the Job Store if restarting is true ($mode)" in {
        ejea = helper.buildEJEA(restarting = true)(RequestingRestartCheckToken)
        ejea ! JobTokenDispensed

        helper.jobStoreProbe.expectMsgPF(max = awaitTimeout, hint = "Awaiting job store lookup") {
          case QueryJobCompletion(jobKey, taskOutputs) =>
            validateJobStoreKey(jobKey)
            taskOutputs should be(helper.call.outputPorts.toSeq)
        }
        helper.bjeaProbe.expectNoMessage(awaitAlmostNothing)
        helper.jobHashingInitializations shouldBe NothingYet
        ejea.stateName should be(CheckingJobStore)
      }


      s"bypass the Job Store and request output store to start preparing the job for running or call caching ($mode)" in {
        ejea = helper.buildEJEA(restarting = false)
        ejea ! JobTokenDispensed

        helper.replyToProbe.expectMsg(max = awaitTimeout, hint = "Awaiting JobStarting message", JobStarting(helper.jobDescriptorKey))
        helper.replyToProbe.expectMsg(max = awaitTimeout, hint = "Awaiting RequestValueStore message", RequestValueStore)
        ejea.stateName should be(WaitingForValueStore)
      }
    }
  }
}
