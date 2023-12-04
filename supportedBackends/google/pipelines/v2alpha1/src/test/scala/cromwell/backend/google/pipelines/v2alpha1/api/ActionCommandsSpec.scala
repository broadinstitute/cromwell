package cromwell.backend.google.pipelines.v2alpha1.api

import java.nio.file.Path

import common.assertion.CromwellTimeoutSpec
import cromwell.backend.google.pipelines.common.PipelinesApiConfigurationAttributes.GcsTransferConfiguration
import cromwell.backend.google.pipelines.common.action.ActionCommands._
import cromwell.filesystems.gcs.GcsPath
import eu.timepit.refined.refineMV
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import common.mock.MockSugar

class ActionCommandsSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers with MockSugar {
  behavior of "ActionCommands"

  it should "inject project flag when request fails because of requester pays" in {
    val path = GcsPath(
      mock[Path],
      mock[com.google.api.services.storage.Storage],
      mock[com.google.cloud.storage.Storage],
      "my-project"
    )
    val recovered = recoverRequesterPaysError(path) { flag =>
      s"flag is $flag"
    }

    recovered shouldBe """flag is  > gsutil_output.txt 2>&1
                         |# Record the exit code of the gsutil command without project flag
                         |RC_GSUTIL=$?
                         |if [ "$RC_GSUTIL" != "0" ]; then
                         |  printf '%s %s\n' "$(date -u '+%Y/%m/%d %H:%M:%S')" flag\ is\ \ failed
                         |  # Print the reason of the failure
                         |  cat gsutil_output.txt
                         |
                         |  # Check if it matches the BucketIsRequesterPaysErrorMessage
                         |  if grep -q "requester pays bucket but no user project" gsutil_output.txt; then
                         |    printf '%s %s\n' "$(date -u '+%Y/%m/%d %H:%M:%S')" Retrying\ with\ user\ project
                         |    flag is -u my-project
                         |  else
                         |    exit "$RC_GSUTIL"
                         |  fi
                         |else
                         |  exit 0
                         |fi""".stripMargin
  }

  it should "use GcsTransferConfiguration to set the number of localization retries" in {
    implicit val gcsTransferConfiguration: GcsTransferConfiguration =
      GcsTransferConfiguration(transferAttempts = refineMV(31380), parallelCompositeUploadThreshold = "0")
    retry(
      "I'm very flaky"
    ) shouldBe """for i in $(seq 31380); do
                 |  (
                 |    I'm very flaky
                 |  )
                 |  RC=$?
                 |  if [ "$RC" = "0" ]; then
                 |    break
                 |  fi
                 |  if [ $i -lt 31380 ]; then
                 |    printf '%s %s\n' "$(date -u '+%Y/%m/%d %H:%M:%S')" Waiting\ 5\ seconds\ and\ retrying
                 |    sleep 5
                 |  fi
                 |done
                 |exit "$RC"""".stripMargin
  }
}
