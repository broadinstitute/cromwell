package cromwell.engine.backend.jes

import java.text.SimpleDateFormat
import java.util.Date

import com.google.api.services.genomics.model.Operation
import scala.collection.JavaConverters._

object OperationMetadata {
  def apply(operation: Operation): OperationMetadata = {
    val metadata = operation.getMetadata.asScala
    new OperationMetadata(metadata("createTime").asInstanceOf[String].toDate,
      metadata.get("startTime").map(_.asInstanceOf[String].toDate),
      metadata.get("endTime").map(_.asInstanceOf[String].toDate))
  }

  val Iso8601DateFormat = new SimpleDateFormat("yyyy-MM-dd'T'HH:mm:ss.SSSXXX")

  implicit class EnhancedString(val string: String) extends AnyVal {
    def toDate: Date = Iso8601DateFormat.parse(string.asInstanceOf[String])
  }
}

case class OperationMetadata(created: Date, started: Option[Date], finished: Option[Date])
