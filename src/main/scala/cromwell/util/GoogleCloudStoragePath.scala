package cromwell.util

import scala.util.{Failure, Success, Try}

case class GoogleCloudStoragePath(bucket: String, objectName: String) {
  override def toString() = {
    "gs://" + bucket + "/" + objectName
  }

  def +(value: String): GoogleCloudStoragePath = {
    if (objectName.endsWith("/") || value.startsWith("/")) {
      GoogleCloudStoragePath(bucket, objectName + value)
    } else {
      GoogleCloudStoragePath(bucket, objectName + "/" + value)
    }
  }
}

object GoogleCloudStoragePath {
  def apply(value: String): GoogleCloudStoragePath = {
    tryParse(value) match {
      case Success(gcsPath) => gcsPath
      case Failure(e) => throw e
    }
  }

  def tryParse(value: String): Try[GoogleCloudStoragePath] = {
    val gsUriRegex = """gs://([^/]*)/(.*)""".r
    value match {
      case gsUriRegex(bucket, objectName) => Success(GoogleCloudStoragePath(bucket, objectName))
      case _ => Failure(new IllegalArgumentException())
    }
  }
}