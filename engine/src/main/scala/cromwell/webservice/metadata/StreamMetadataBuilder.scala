package cromwell.webservice.metadata

import cats.effect.IO
import com.typesafe.scalalogging.StrictLogging
import cromwell.services.MetadataServicesStore
import cromwell.services.metadata.impl.MetadataDatabaseAccess
import cromwell.services.metadata.{MetadataEvent, MetadataKey, MetadataValue, _}
import cromwell.webservice.metadata.MetadataComponent._
import fs2.interop.reactivestreams._
import spray.json._

import scala.concurrent.ExecutionContext

object StreamMetadataBuilder extends StrictLogging with MetadataDatabaseAccess with MetadataServicesStore {
  
  def build(query: MetadataQuery)(implicit ec: ExecutionContext): IO[JsObject] = {

    val eventStream = queryToEventStream(query)
    
    // Workflow Id is not pushed as a metadata event, so manually add it here
    val workflowIdEvent = MetadataEvent(MetadataKey(query.workflowId, None, "id"), MetadataValue(query.workflowId.toString))

    (eventStream ++ fs2.Stream.emit(workflowIdEvent))
      .foldMap(toMetadataComponent(query, ec))
      .compile
      .toVector
      .map(_.head)
      .map(_.toJson.asJsObject)
  }

  private def queryToEventStream(query: MetadataQuery)(implicit ec: ExecutionContext) = {
    implicit val contextShift = IO.contextShift(ec)
    val subscriber = StreamSubscriber[IO, MetadataEvent].unsafeRunSync()
    val databasePublisher = streamedQueryMetadataEvents(query)
    databasePublisher.subscribe(subscriber)
    subscriber.stream
  }
}
