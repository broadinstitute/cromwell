package cromwell.webservice.metadata

import cats.effect.IO
import cromwell.services.metadata.{MetadataEvent, MetadataKey, MetadataValue}
import slick.basic.DatabasePublisher
import fs2.interop.reactivestreams._
import MetadataComponent._
import cromwell.core.WorkflowId

import scala.concurrent.ExecutionContext
import spray.json._

object StreamMetadataBuilder {
  def build(databasePublisher: DatabasePublisher[MetadataEvent], workflowId: WorkflowId)(implicit ec: ExecutionContext): JsObject = {
    implicit val contextShift = IO.contextShift(ec)
    val subscriber = StreamSubscriber[IO, MetadataEvent].unsafeRunSync()
    databasePublisher.subscribe(subscriber)
    
    val workflowIdEvent = MetadataEvent(
      MetadataKey(workflowId, None, "id"),
      MetadataValue(workflowId.toString)
    )
    
    val stream = (subscriber.stream ++ fs2.Stream.emit(workflowIdEvent))
      .foldMap(toMetadataComponent(Map.empty))
      .map(_.toJson.asJsObject)

    stream.compile.toVector.unsafeRunSync().head
  }
}
