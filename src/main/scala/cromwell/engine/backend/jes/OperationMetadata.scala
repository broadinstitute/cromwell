package cromwell.engine.backend.jes

import com.typesafe.scalalogging.LazyLogging
import com.google.api.services.genomics.model.Operation
import org.joda.time.DateTime
import scala.collection.JavaConverters._

object OperationMetadata extends LazyLogging {
  def apply(operation: Operation): OperationMetadata = {
    val metadata = operation.getMetadata.asScala
    new OperationMetadata(metadata("createTime").asInstanceOf[String].toDate,
      metadata.get("startTime").map(_.asInstanceOf[String].toDate),
      metadata.get("endTime").map(_.asInstanceOf[String].toDate))
  }

  val defaultDate = DateTime.parse("0")

  implicit class EnhancedString(val string: String) extends AnyVal {
    def toDate: DateTime = try {

      DateTime.parse(string.asInstanceOf[String])
    } catch {
      case e: Exception =>
        logger.error(s"Unexpected date string: '$string'.")
        defaultDate
    }
  }
}

case class OperationMetadata(created: DateTime, started: Option[DateTime], finished: Option[DateTime])
