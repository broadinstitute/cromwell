package cromwell.engine.backend.jes

import com.google.api.services.genomics.model.Operation
import scala.collection.JavaConverters._

object OperationMetadata {
  def apply(operation: Operation): OperationMetadata = {
    val metadata = operation.getMetadata.asScala
    new OperationMetadata(metadata("createTime").asInstanceOf[String],
      metadata.get("startTime").map {_.asInstanceOf[String]},
      metadata.get("endTime").map {_.asInstanceOf[String]})
  }
}

// FIXME: These dates aren't strings
// FIXME: What's up w/ the options? Clarify w/ the ADT as well
case class OperationMetadata(created: String, started: Option[String], finished: Option[String])
