package cromwell.backend.impl.bcs

import org.scalatest.TryValues._

class RunStatusFactorySpec extends BcsTestUtilSpec {
  behavior of s"RunStatusFactorySpec"
  import RunStatus._

  it should "get a right status from factory" in {
    RunStatusFactory.getStatus(jobId, "Waiting").success.value shouldBe a [Waiting]
    RunStatusFactory.getStatus(jobId, "Running").success.value shouldBe a [Running]
    RunStatusFactory.getStatus(jobId, "Finished").success.value shouldBe a [Finished]
    RunStatusFactory.getStatus(jobId, "Failed").success.value shouldBe a [Failed]
    RunStatusFactory.getStatus(jobId, "Stopped").success.value shouldBe a [Stopped]
  }

  it should "" in {

  }
}
