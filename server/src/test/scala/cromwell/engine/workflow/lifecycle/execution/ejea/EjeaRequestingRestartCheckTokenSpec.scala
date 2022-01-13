package cromwell.engine.workflow.lifecycle.execution.ejea

import cromwell.engine.workflow.lifecycle.execution.job.EngineJobExecutionActor._
import cromwell.engine.workflow.tokens.JobTokenDispenserActor.JobTokenDispensed
import cromwell.jobstore.JobStoreActor.QueryJobCompletion
import org.scalatest.concurrent.Eventually

class EjeaRequestingRestartCheckTokenSpec extends EngineJobExecutionActorSpec with CanValidateJobStoreKey with Eventually {

  override implicit val stateUnderTest: EngineJobExecutionActorState = RequestingRestartCheckToken

  "An EJEA in the RequestingRestartTokenCheck state" should {

    s"do nothing when denied a token (with restarting=true)" in {
      ejea = helper.buildEJEA(restarting = true)
      helper.jobRestartCheckTokenDispenserProbe.expectNoMessage(max = awaitAlmostNothing)
      ejea.stateName should be(RequestingRestartCheckToken)
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
    }
  }
}
