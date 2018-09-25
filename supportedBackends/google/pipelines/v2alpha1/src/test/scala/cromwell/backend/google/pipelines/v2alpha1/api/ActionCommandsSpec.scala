package cromwell.backend.google.pipelines.v2alpha1.api

import java.nio.file.Path

import cromwell.backend.google.pipelines.common.PipelinesApiAttributes.LocalizationConfiguration
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
    
    recovered shouldBe "flag is  2> gsutil_output.txt; RC_GSUTIL=$?; if [ \"$RC_GSUTIL\" = \"1\" ]; then\n grep \"Bucket is requester pays bucket but no user project provided.\" gsutil_output.txt && echo \"Retrying with user project\"; flag is -u my-project; fi "
  }
  
  it should "use LocalizationConfiguration to set the number of localization retries" in {
    implicit val localizationConfiguration = LocalizationConfiguration(refineMV(31380))
    retry("I'm very flaky") shouldBe "retry() { for i in `seq 31380`; do I'm very flaky; RC=$?; if [ \"$RC\" = \"0\" ]; then break; fi; sleep 5; done; return \"$RC\"; }; retry"
  }
}
