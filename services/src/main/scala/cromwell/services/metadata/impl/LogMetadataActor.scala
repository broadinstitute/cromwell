package cromwell.services.metadata.impl

import akka.actor.Actor
import cromwell.services.metadata.MetadataEvent
import cromwell.services.metadata.MetadataService.PutMetadataAction
import io.circe.Printer

case class LogMetadataActor(logger: String => Unit) extends Actor {
  import LogMetadataActor._
  override def receive: Receive = {
    case PutMetadataAction(events: Iterable[MetadataEvent]) => events foreach (logger compose metadataToString)
  }
}

object LogMetadataActor {
  import io.circe.generic.auto._
  import io.circe.syntax._

  final val compact: Printer = Printer(
    preserveOrder = true,
    dropNullValues = true,
    indent = ""
  )

  def metadataToString: MetadataEvent => String = e =>
    s"${e.asJson.pretty(compact)}"
}
