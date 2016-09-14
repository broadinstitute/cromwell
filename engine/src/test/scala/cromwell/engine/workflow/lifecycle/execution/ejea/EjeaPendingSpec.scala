package cromwell.engine.workflow.lifecycle.execution.ejea

import cromwell.engine.workflow.lifecycle.execution.EngineJobExecutionActor.{CheckingJobStore, EngineJobExecutionActorState, Execute, Pending, PreparingJob}
import cromwell.engine.workflow.lifecycle.execution.JobPreparationActor
import cromwell.jobstore.JobStoreActor.QueryJobCompletion
import org.scalatest.concurrent.Eventually

class EjeaPendingSpec extends EngineJobExecutionActorSpec with CanValidateJobStoreKey with Eventually {

  override implicit val stateUnderTest: EngineJobExecutionActorState = Pending

  "An EJEA in the Pending state" should {

    CallCachingModes foreach { mode =>
      s"check against the Job Store if restarting is true ($mode)" in {
        ejea = helper.buildEJEA(restarting = true)
        ejea ! Execute

        helper.jobStoreProbe.expectMsgPF(max = awaitTimeout, hint = "Awaiting job store lookup") {
          case QueryJobCompletion(jobKey, taskOutputs) =>
            validateJobStoreKey(jobKey)
            taskOutputs should be(helper.task.outputs)
        }
        helper.bjeaProbe.expectNoMsg(awaitAlmostNothing)
        helper.jobHashingInitializations shouldBe NothingYet
        ejea.stateName should be(CheckingJobStore)
      }


      s"bypass the Job Store and start preparing the job for running or call caching ($mode)" in {
        ejea = helper.buildEJEA(restarting = false)
        ejea ! Execute

        helper.jobPreparationProbe.expectMsg(max = awaitTimeout, hint = "Awaiting job preparation", JobPreparationActor.Start)
        helper.jobStoreProbe.expectNoMsg(awaitAlmostNothing)
        ejea.stateName should be(PreparingJob)
      }
    }
  }
}
