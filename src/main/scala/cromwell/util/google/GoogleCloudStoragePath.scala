package cromwell.util.google

import scala.util.{Failure, Success, Try}

case class GoogleCloudStoragePath(bucket: String, objectName: String) {
  override def toString = {
    "gs://" + bucket + "/" + objectName
  }
}

object GoogleCloudStoragePath {
  def apply(value: String): GoogleCloudStoragePath = {
    parse(value) match {
      case Success(gcsPath) => gcsPath
      case Failure(e) => throw e
    }
  }

  def parse(value: String): Try[GoogleCloudStoragePath] = {
    val gsUriRegex = """gs://([^/]*)/(.*)""".r
    value match {
      case gsUriRegex(bucket, objectName) => Success(GoogleCloudStoragePath(bucket, objectName))
      case _ => Failure(new IllegalArgumentException())
    }
  }
}