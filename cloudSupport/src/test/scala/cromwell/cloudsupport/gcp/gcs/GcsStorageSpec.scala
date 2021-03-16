package cromwell.cloudsupport.gcp.gcs

import com.google.api.client.googleapis.batch.BatchRequest
import com.google.api.client.googleapis.batch.json.JsonBatchCallback
import com.google.api.client.googleapis.json.GoogleJsonError
import com.google.api.client.http.{HttpHeaders, HttpRequest, HttpRequestInitializer}
import com.google.api.client.json.jackson2.JacksonFactory
import com.google.api.gax.retrying.RetrySettings
import com.google.api.services.storage.{Storage, StorageRequest}
import com.google.api.services.storage.model.RewriteResponse
import com.google.auth.oauth2.GoogleCredentials
import com.google.cloud.NoCredentials
import common.assertion.CromwellTimeoutSpec
import cromwell.cloudsupport.gcp.GoogleConfiguration
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers


class GcsStorageSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers {

  behavior of "GcsStorage"

  it should "build the default cloud storage configuration" ignore {
    val configuration = GcsStorage.DefaultCloudStorageConfiguration
    configuration.permitEmptyPathComponents should be(true)
    configuration.stripPrefixSlash should be(true)
    configuration.usePseudoDirectories should be(true)
  }

  it should "build gcs storage" ignore {
    val configuration = GcsStorage.gcsStorage("gcs-storage-spec", NoCredentials.getInstance(), RetrySettings.newBuilder().build())
    configuration.getApplicationName should be("gcs-storage-spec")
  }

  it should "build gcs storage options" ignore {
    val projectOption = Option("custom-project")
    val retrySettings = RetrySettings.newBuilder().setMaxAttempts(2).build()
    val storageOptions = GcsStorage.gcsStorageOptions(NoCredentials.getInstance(), retrySettings, projectOption)
    storageOptions.getProjectId should be("custom-project")
    storageOptions.getRetrySettings.getMaxAttempts should be(2)
  }

  private val httpRequestInitializer = new HttpRequestInitializer {
    override def initialize(request: HttpRequest): Unit = {
      request.setConnectTimeout(GoogleConfiguration.DefaultConnectionTimeout.toMillis.toInt)
      request.setReadTimeout(GoogleConfiguration.DefaultReadTimeout.toMillis.toInt)
      ()
    }
  }

  it should "handle large GCS batches" in {
    // Inside the batch, per entry, embed the credentials
    val apiStorage: Storage =
      GcsStorage.gcsStorage(
        "gcs-storage-spec",
        GoogleCredentials.getApplicationDefault(),
        RetrySettings.newBuilder().setMaxAttempts(1).build(),
      )
    // For the outer batch, do not set the credentials
    val batchStorage: Storage = new Storage.Builder(
      GcsStorage.HttpTransport,
      JacksonFactory.getDefaultInstance,
      httpRequestInitializer
    ).setApplicationName("gcs-storage-spec").build()

    val batchRequest: BatchRequest = batchStorage.batch()

    val bucket: String = "kshakir-large-num-files-ok2delete"
    val blobDir: String = "0" * 1//210
    for {
      i <- 0 until 3001
    } {
      val blobName: String = f"$blobDir%s/$i%04d.txt"

      val rewrite: StorageRequest[RewriteResponse] = apiStorage.objects()
        .rewrite(bucket, "empty.txt", bucket, s"src/$blobName", null)
        //.rewrite(bucket, s"src/$blobName", bucket, s"dst/$blobName", null)

      rewrite.queue(batchRequest, new JsonBatchCallback[RewriteResponse] {
        override def onFailure(e: GoogleJsonError, responseHeaders: HttpHeaders): Unit = {
          println(s"Failure for $i:\n$e")
        }

        override def onSuccess(t: RewriteResponse, responseHeaders: HttpHeaders): Unit = {
          println(s"Success for $i")
        }
      })
    }

    batchRequest.execute()
  }

}
