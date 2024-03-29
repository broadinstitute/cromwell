package cromwell.backend.impl.sfs.config

import common.assertion.CromwellTimeoutSpec
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

class ConfigAsyncJobExecutionActorSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers {

  behavior of "getJobRegex"

  val successfulTestCases = Map(
    ("Job <03957>... blah blah blah", "Job <(\\d+)>.*") -> "03957",
    ("""mxq_group_id=...
       |mxq_group_name=...
       |mxq_job_id=16988030
       |""".stripMargin,
     "mxq_job_id=(\\d+)"
    ) -> "16988030"
  )

  successfulTestCases foreach { case ((fileContent, jobIdRegex), expectedJobId) =>
    it should s"find the right job ID $expectedJobId using the regex $jobIdRegex" in {
      DispatchedConfigAsyncJobExecutionActor.getJob(fileContent, null, jobIdRegex).jobId should be(expectedJobId)
    }
  }
}
