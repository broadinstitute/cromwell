package cromwell.binding.values

import cromwell.util.GoogleCloudStoragePath
import scala.util.{Success, Try}

object WdlGcsObject {
  def apply(value: String): WdlGcsObject = WdlGcsObject(GoogleCloudStoragePath(value))
  def apply(bucket: String, objectName: String): WdlGcsObject = WdlGcsObject(GoogleCloudStoragePath(bucket, objectName))
}

/**
 * Equivalent to WdlFile, for GCS objects.
 * @param value Path for the GCS object.
 */
case class WdlGcsObject(value: GoogleCloudStoragePath) extends WdlFileLike {

  /**
   * Addition undefined for GCS paths
   * @param rhs Ignored
   * @return Failure.
   */
  override def add(rhs: WdlValue): Try[WdlValue] = invalid(s"$value + $rhs")

  override def equals(rhs: WdlValue): Try[WdlBoolean] = {
    rhs match {
      case r: WdlGcsObject => Success(WdlBoolean(value.equals(r.value)))
      case _ => invalid(s"$value == $rhs")
    }
  }

  override def asString = value.toString
}
