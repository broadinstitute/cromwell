package cromwell.util

import java.net.URI

case class GoogleCloudStoragePath(bucket: String, objectName: String)

object GoogleCloudStoragePath {
  def apply(value: String): GoogleCloudStoragePath = {
    val gsUriRegex = """gs://([^/]*)/(.*)""".r
    value match {
      case gsUriRegex(bucket, objectName) => GoogleCloudStoragePath(bucket, objectName)
      case _ => throw new IllegalArgumentException()
    }
  }
}