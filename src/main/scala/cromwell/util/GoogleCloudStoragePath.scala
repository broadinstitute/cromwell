package cromwell.util

import scala.util.{Failure, Success, Try}

case class GoogleCloudStoragePath(bucket: String, objectName: String) {
  override def toString() = {
    "gs://" + bucket + "/" + objectName
  }

  def +(value: String): GoogleCloudStoragePath = {
    val trimmedObjectName = if (objectName.endsWith("/")) { objectName.substring(0, objectName.length - 1)} else { objectName }
    val trimmedValue = if (value.startsWith("/")) { value.substring(1, value.length)} else { value }

    GoogleCloudStoragePath(s"gs://$bucket/$trimmedObjectName/$trimmedValue")
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
      case _ => Failure(new IllegalArgumentException(value))
    }
  }
}