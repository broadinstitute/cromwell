package cromwell.filesystems.gcs

import java.util.concurrent.Executors

import com.google.api.gax.retrying.RetrySettings
import com.google.cloud.NoCredentials
import com.google.cloud.storage.contrib.nio.CloudStorageConfiguration
import cromwell.cloudsupport.gcp.gcs.GcsStorage
import cromwell.core.WorkflowOptions

import scala.concurrent.ExecutionContext

object MockGcsPathBuilder {
  implicit val ec = ExecutionContext.fromExecutor(Executors.newSingleThreadExecutor())
  private def makeStorageOptions(project: Option[String] = Option("cromwell-test")) =
    GcsStorage.gcsStorageOptions(NoCredentials.getInstance(), RetrySettings.newBuilder().build(), project)
  private val storageOptions = makeStorageOptions()
  private val apiStorage = GcsStorage.gcsStorage("cromwell-test-app", storageOptions)

  lazy val instance = new GcsPathBuilder(apiStorage, CloudStorageConfiguration.DEFAULT, storageOptions)

  def withOptions(workflowOptions: WorkflowOptions) = {
    val customStorageOptions = makeStorageOptions(workflowOptions.get("google_project").toOption)
    new GcsPathBuilder(apiStorage, GcsStorage.DefaultCloudStorageConfiguration, customStorageOptions)
  }
}
