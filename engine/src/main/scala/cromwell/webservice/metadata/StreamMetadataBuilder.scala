package cromwell.webservice.metadata

import cats.effect.IO
import cats.effect.IO._
import cats.instances.vector._
import cats.instances.tuple._
import cats.syntax.parallel._
import common.validation.Validation._
import com.typesafe.scalalogging.LazyLogging
import cromwell.core.{WorkflowId, WorkflowMetadataKeys}
import cromwell.services.MetadataServicesStore
import cromwell.services.metadata._
import cromwell.services.metadata.impl.StreamMetadataDatabaseAccess
import cromwell.webservice.metadata.MetadataComponent._
import fs2.interop.reactivestreams._
import fs2.{Pipe, Stream => FStream}
import org.reactivestreams.Publisher
import spray.json._

import scala.concurrent.ExecutionContext
import scala.concurrent.duration._

class StreamMetadataBuilder(requestTimeout: FiniteDuration) extends StreamMetadataDatabaseAccess with MetadataServicesStore with LazyLogging {
  def buildQuery(query: MetadataQuery,
                 eventToComponent: MetadataEvent => MetadataComponent,
                 extraComponents: Seq[MetadataComponent] = Seq.empty)
                (implicit ec: ExecutionContext): IO[JsObject] = {
    implicit val contextShift = cats.effect.IO.contextShift(ec)

    // Given a vector of workflow ids, build the corresponding jsons in parallel using the same query
    def buildSubWorkflows(ids: Vector[WorkflowId]): IO[Map[WorkflowId, JsObject]] = {
      ids.parTraverse[IO, IO.Par, (WorkflowId, JsObject)](id =>
        buildQuery(query.copy(workflowId = id), eventToComponent, extraComponents).map(id -> _)
      ).map(_.toMap)
    }
    
    // Transform a MetadataEvent to a MetadataComponent, and optionally add it to the sub workflow id list if necessary
    def mapEvent(subWorkflowIds: Vector[WorkflowId], event: MetadataEvent): (Vector[WorkflowId], MetadataComponent) = {
      if (event.isSubWorkflowId)
        (event.extractSubWorkflowId.map(subWorkflowIds.:+(_)).getOrElse(subWorkflowIds), eventToComponent(event))
      else
        (subWorkflowIds, eventToComponent(event))
    }
    
    // Run the stream by mapping metadata events to metadata components while accumulating sub workflow ids. Then fold everything together
    def runStream(stream: FStream[IO, MetadataEvent]): IO[(Vector[WorkflowId], MetadataComponent)] = stream
      .mapAccumulate[Vector[WorkflowId], MetadataComponent](Vector.empty[WorkflowId])(mapEvent)
      .append(fs2.Stream.emits(extraComponents.map(Vector.empty[WorkflowId] -> _)))
      .compile
      .foldMonoid

    for {
      stream <- queryToEventStream(query)
      subWorkflowIdsAndComponent <- runStream(stream)
      (subWorkflowIds, component) = subWorkflowIdsAndComponent
      subWorkflowJsons <- buildSubWorkflows(subWorkflowIds)
      componentWriter = new MetadataComponentJsonWriter(subWorkflowJsons)
    } yield component.toJson(componentWriter).asJsObject
  }

  def workflowIdComponent(id: WorkflowId): MetadataComponent = MetadataObject(Map("id" -> MetadataPrimitive(MetadataValue(id.toString))))
  lazy val workflowOutputsComponent: MetadataComponent = MetadataObject(Map(WorkflowMetadataKeys.Outputs -> MetadataEmptyComponent))

  def queryToEventStream(query: MetadataQuery)(implicit ec: ExecutionContext): IO[fs2.Stream[IO, MetadataEvent]] = {
    for {
      publisher <- queryMetadataEvents(query).toIO("Invalid metadata query")
      stream <- publisherToEventStream(publisher, Option(query.workflowId))
    } yield stream
  }

  /**
    * Builds the json corresponding to the provided MetadataQuery
    */
  def workflowMetadataQuery(query: MetadataQuery)(implicit ec: ExecutionContext): IO[JsObject] = {
    buildQuery(query, toMetadataComponent, List(workflowIdComponent(query.workflowId)))
  }

  /**
    * Creates an fs2 stream from a Publisher of MetadataEvents. The requestTimeout is applied, meaning the stream will be 
    * terminated if processing goes over requestTimeout0
    * @param eventsPublisher publisher of metadata events
    * @param workflowId optional workflow id (for logging in case of timetout)
    */
  private def publisherToEventStream(eventsPublisher: Publisher[MetadataEvent], workflowId: Option[WorkflowId])(implicit ec: ExecutionContext): IO[FStream[IO, MetadataEvent]] = {
    implicit val timer = cats.effect.IO.timer(ec)
    implicit val contextShift = IO.contextShift(ec)

    // Instantiating a subscriber is wrapped in an IO so unwrap it here, connect it to the publisher and get the stream
    StreamSubscriber[IO, MetadataEvent] map { subscriber =>
      val databasePublisher = eventsPublisher
      databasePublisher.subscribe(subscriber)

      subscriber
        .stream
        .interruptWhen[IO]((FStream.sleep_[IO](requestTimeout) ++ FStream(true)) through logTimeout(workflowId))
    }
  }

  /**
    * Utility function to log when a request is forcefully timed out
    */
  private def logTimeout[A](workflowId : Option[WorkflowId]) : Pipe[IO,A,A] = { in =>
    in.evalMap({ element =>
      IO {
        logger.warn(s"Timed out building metadata for workflow: ${workflowId.getOrElse("Unknown")}")
        element
      }
    })
  }
}
