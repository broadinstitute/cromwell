package cromwell.backend.impl.bcs

import org.scalatest.TryValues._



final class RunStatusFactorySpec extends BcsTestUtilSpec {
  behavior of s"RunStatusFactorySpec"

  private case class Status(str: String,
                          isRunningOrComplete: Boolean,
                          terminated: Boolean)

  strToClasses foreach { status =>
    it should behave like verifyStatus(status)
  }

  def verifyStatus(status: Status) = {
    it should s"have correct status: ${status.str}" in withClue(status.str) {
      val s = RunStatusFactory.getStatus(jobId, status.str).success.value
      s.status shouldEqual status.str
      s.isRunningOrComplete shouldEqual status.isRunningOrComplete
      s.isTerminated shouldEqual status.terminated
    }
  }

  private def strToClasses = Seq(
    Status("Waiting", false, false),
    Status("Running", true, false),
    Status("Finished", true, true),
    Status("Stopped", true, true),
    Status("Failed", true, true)
  )
}
