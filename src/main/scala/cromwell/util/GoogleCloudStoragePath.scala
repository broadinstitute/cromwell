package cromwell.util

import java.net.URI

object GoogleCloudStoragePath {
  def apply(value: String): GoogleCloudStoragePath = {
    val gsUriRegex = """gs://([^/]*)/(.*)""".r
    value match {
      case gsUriRegex(bucket, objectName) => GoogleCloudStoragePath(bucket, objectName)
      // TODO: error handle:
      case _ => throw new IllegalArgumentException()
    }
  }

  def formUriString(bucket: String, objectName: String) = "gs://" + bucket + "/" + objectName

  def asURI(gcsUri: GoogleCloudStoragePath) = new URI(formUriString(gcsUri.bucket, gcsUri.objectName))

  def asURI(bucket: String, objectName: String) = new URI(formUriString(bucket, objectName))
}

case class GoogleCloudStoragePath(bucket: String, objectName: String)