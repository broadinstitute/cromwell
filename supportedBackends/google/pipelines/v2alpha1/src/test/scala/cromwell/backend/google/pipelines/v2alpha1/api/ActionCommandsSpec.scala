package cromwell.backend.google.pipelines.v2alpha1.api

import java.nio.file.Path

import cromwell.backend.google.pipelines.common.PipelinesApiConfigurationAttributes.LocalizationConfiguration
import cromwell.filesystems.gcs.GcsPath
import org.scalatest.{FlatSpec, Matchers}
import org.specs2.mock.Mockito
import eu.timepit.refined.refineMV
import ActionCommands._

class ActionCommandsSpec extends FlatSpec with Matchers with Mockito {
  behavior of "ActionCommands"
  
  it should "inject project flag when request fails because of requester pays" in {
    val path = GcsPath(any[Path], any[com.google.api.services.storage.Storage], any[com.google.cloud.storage.Storage], "my-project")
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
                         |  if grep -q "Bucket is requester pays bucket but no user project provided." gsutil_output.txt; then
                         |    printf '%s %s\n' "$(date -u '+%Y/%m/%d %H:%M:%S')" Retrying\ with\ user\ project
                         |    flag is -u my-project
                         |  else
                         |    exit "$RC_GSUTIL"
                         |  fi
                         |else
                         |  exit 0
                         |fi""".stripMargin
  }
  
  it should "use LocalizationConfiguration to set the number of localization retries" in {
    implicit val localizationConfiguration = LocalizationConfiguration(refineMV(31380))
    retry("I'm very flaky") shouldBe """for i in $(seq 31380); do
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
