package cromwell.util.google

import scala.util.{Failure, Success, Try}
import scala.language.implicitConversions

/** Represents a Google Cloud Storage path, like gs://bucket/path/to/object.txt
  *
  * @param bucket - should adhere to https://cloud.google.com/storage/docs/bucket-naming?hl=en#requirements
  * @param objectName
  */
case class GcsPath(bucket: String, objectName: String) {
  override def toString = {
    "gs://" + bucket + "/" + objectName
  }
}

object GcsPath {

  implicit def toGcsPath(str: String): GcsPath = GcsPath(str)

  def apply(value: String): GcsPath = {
    parse(value) match {
      case Success(gcsPath) => gcsPath
      case Failure(e) => throw e
    }
  }

  def parse(value: String): Try[GcsPath] = {
    val gsUriRegex = """gs://([^/]*)(.*)""".r
    value match {
      case gsUriRegex(bucket, objectName) => Success(GcsPath(bucket, objectName.stripPrefix("/")))
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