package cromwell.cloudsupport.gcp.gcs

import com.google.api.gax.retrying.RetrySettings
import com.google.cloud.NoCredentials
import common.assertion.CromwellTimeoutSpec
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers


class GcsStorageSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers {

  behavior of "GcsStorage"

  it should "build the default cloud storage configuration" in {
    val configuration = GcsStorage.DefaultCloudStorageConfiguration
    configuration.permitEmptyPathComponents should be(true)
    configuration.stripPrefixSlash should be(true)
    configuration.usePseudoDirectories should be(true)
  }

  it should "build gcs storage" in {
    val configuration = GcsStorage.gcsStorage("gcs-storage-spec", NoCredentials.getInstance(), RetrySettings.newBuilder().build())
    configuration.getApplicationName should be("gcs-storage-spec")
  }

  it should "build gcs storage options" in {
    val projectOption = Option("custom-project")
    val retrySettings = RetrySettings.newBuilder().setMaxAttempts(2).build()
    val storageOptions = GcsStorage.gcsStorageOptions(NoCredentials.getInstance(), retrySettings, projectOption)
    storageOptions.getProjectId should be("custom-project")
    storageOptions.getRetrySettings.getMaxAttempts should be(2)
  }

}
