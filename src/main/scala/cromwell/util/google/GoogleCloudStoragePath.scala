package cromwell.util.google

import scala.util.{Failure, Success, Try}
import scala.language.implicitConversions

case class GoogleCloudStoragePath(bucket: String, objectName: String) {
  override def toString = {
    "gs://" + bucket + "/" + objectName
  }
}

object GoogleCloudStoragePath {

  implicit def toGcsPath(str: String): GoogleCloudStoragePath = GoogleCloudStoragePath(str)

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
      case _ => Failure(new IllegalArgumentException(s"Not a valid Google Cloud Storage URI: $value"))
    }
  }

  def parseBucket(value: String): Try[String] = {
    val gsUriRegex = """gs://([^/]*)""".r
    value match {
      case gsUriRegex(bucket) => Success(bucket)
      case _ => Failure(new IllegalArgumentException(s"Not a valid Google Cloud Storage URI: $value"))
    }
  }
}