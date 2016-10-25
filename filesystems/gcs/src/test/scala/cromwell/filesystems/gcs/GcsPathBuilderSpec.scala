package cromwell.filesystems.gcs

import com.google.cloud.RetryParams
import com.google.cloud.storage.contrib.nio.CloudStorageConfiguration
import cromwell.core.path.CustomRetryParams
import cromwell.core.path.proxy.RetryableFileSystemProviderProxy
import cromwell.core.{TestKitSuite, WorkflowOptions}
import cromwell.filesystems.gcs.auth.GoogleAuthMode
import org.scalatest.{FlatSpecLike, Matchers}

class GcsPathBuilderSpec extends TestKitSuite with FlatSpecLike with Matchers {

  implicit val as = system

  behavior of "GcsPathBuilderSpec"

  it should "create a path with a retryable provider" in {
    val retryablePathBuilder = new RetryableGcsPathBuilder(
      GoogleAuthMode.NoAuthMode,
      RetryParams.defaultInstance(),
      CustomRetryParams.Default,
      CloudStorageConfiguration.DEFAULT,
      WorkflowOptions.empty
    )

    val path = retryablePathBuilder.build("gs://bucket/object")
    path.isSuccess shouldBe true
    path.get.getFileSystem.provider() shouldBe a[RetryableFileSystemProviderProxy[_]]
  }

}
