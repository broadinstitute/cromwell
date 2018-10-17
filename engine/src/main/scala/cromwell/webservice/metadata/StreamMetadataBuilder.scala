package cromwell.webservice.metadata

import cats.effect.IO
import cats.effect.IO._
import cats.instances.vector._
import cats.instances.tuple._
import cats.syntax.parallel._
import com.typesafe.scalalogging.LazyLogging
import cromwell.core.{WorkflowId, WorkflowMetadataKeys}
import cromwell.services.MetadataServicesStore
import cromwell.services.metadata._
import cromwell.services.metadata.impl.MetadataDatabaseAccess
import cromwell.webservice.metadata.MetadataComponent._
import fs2.interop.reactivestreams._
import fs2.{Pipe, Stream}
import org.reactivestreams.Publisher
import spray.json._

import scala.concurrent.ExecutionContext
import scala.concurrent.duration._

object StreamMetadataBuilder extends MetadataDatabaseAccess with MetadataServicesStore with LazyLogging {
  private def logTimeout[A](workflowId : Option[WorkflowId]) : Pipe[IO,A,A] = { in =>
    in.evalMap({ element => 
      IO { 
        logger.warn(s"Timed out building metadata for workflow: ${workflowId.getOrElse("Unknown")}")
        element
      }
    })
  }

  val DefaultTimeout = 30.hours

  def buildFromPublisher(eventsPublisher: Publisher[MetadataEvent], workflowId: Option[WorkflowId] = None, extraComponents: Seq[MetadataComponent] = Seq.empty)(implicit ec: ExecutionContext): IO[MetadataComponent] = {
    build(publisherToEventStream(eventsPublisher), toMetadataComponent, workflowId, extraComponents)
  }

  def build(eventStream: fs2.Stream[IO, MetadataEvent],
            eventToComponent: MetadataEvent => MetadataComponent,
            workflowId: Option[WorkflowId],
            extraComponents: Seq[MetadataComponent] = Seq.empty)(implicit ec: ExecutionContext): IO[MetadataComponent] = {
    implicit val timer = cats.effect.IO.timer(ec)
    implicit val contextShift = cats.effect.IO.contextShift(ec)

    eventStream
      .interruptWhen[IO]((Stream.sleep_[IO](DefaultTimeout) ++ Stream(true)) through logTimeout(workflowId))
      .map(eventToComponent)
      .append(fs2.Stream.emits(extraComponents))
      .compile
      .foldMonoid
  }
  
  def buildQuery(query: MetadataQuery,
                 eventToComponent: MetadataEvent => MetadataComponent,
                 workflowId: Option[WorkflowId],
                 extraComponents: Seq[MetadataComponent] = Seq.empty)(implicit ec: ExecutionContext): IO[JsObject] = {
    implicit val timer = cats.effect.IO.timer(ec)
    implicit val contextShift = cats.effect.IO.contextShift(ec)

    // Given a vector of workflow ids, build the corresponding jsons in parallel using the same query
    def buildSubWorkflows(ids: Vector[WorkflowId]): IO[Map[WorkflowId, JsObject]] = {
      ids.parTraverse[IO, IO.Par, (WorkflowId, JsObject)](id =>
        buildQuery(query.copy(workflowId = id), eventToComponent, Option(id), extraComponents).map(id -> _)
      ).map(_.toMap)
    }
    
    // Transform a MetadataEvent to a MetadataComponent, and optionally add it to the sub workflow id list if necessary
    def mapEvent(subWorkflowIds: Vector[WorkflowId], event: MetadataEvent): (Vector[WorkflowId], MetadataComponent) = {
      if (event.isSubWorkflowId)
        (event.extractSubWorkflowId.map(subWorkflowIds.:+(_)).getOrElse(subWorkflowIds), eventToComponent(event))
      else
        (subWorkflowIds, eventToComponent(event))
    }

    val streamResult = queryToEventStream(query)
      .interruptWhen[IO]((Stream.sleep_[IO](DefaultTimeout) ++ Stream(true)) through logTimeout(workflowId))
      .mapAccumulate[Vector[WorkflowId], MetadataComponent](Vector.empty[WorkflowId])(mapEvent)
      .append(fs2.Stream.emits(extraComponents.map(Vector.empty[WorkflowId] -> _)))
      .compile
      .foldMonoid

    for {
      subWorkflowIdsAndComponent <- streamResult
      (subWorkflowIds, component) = subWorkflowIdsAndComponent
      subWorkflowJsons <- buildSubWorkflows(subWorkflowIds)
      componentWriter = new MetadataComponentJsonWriter(subWorkflowJsons)
    } yield component.toJson(componentWriter).asJsObject
  }

  def workflowIdComponent(id: WorkflowId): MetadataComponent = MetadataObject(Map("id" -> MetadataPrimitive(MetadataValue(id.toString))))
  lazy val workflowOutputsComponent: MetadataComponent = MetadataObject(Map(WorkflowMetadataKeys.Outputs -> MetadataEmptyComponent))

  def queryToEventStream(query: MetadataQuery)(implicit ec: ExecutionContext): fs2.Stream[IO, MetadataEvent] = {
    publisherToEventStream(queryMetadataEvents(query).valueOr(_ => throw new Exception("blah")))
  }

  def publisherToEventStream(eventsPublisher: Publisher[MetadataEvent])(implicit ec: ExecutionContext): fs2.Stream[IO, MetadataEvent] = {
    implicit val contextShift = IO.contextShift(ec)
    val subscriber = StreamSubscriber[IO, MetadataEvent].unsafeRunSync()
    val databasePublisher = eventsPublisher
    databasePublisher.subscribe(subscriber)
    subscriber.stream
  }

  def workflowMetadataQuery(query: MetadataQuery)(implicit ec: ExecutionContext): IO[JsObject] = {
    buildQuery(query, toMetadataComponent(query, ec), Option(query.workflowId), List(workflowIdComponent(query.workflowId)))
  }

  def workflowOutputs(id: WorkflowId)(implicit ec: ExecutionContext): IO[MetadataComponent] = {
    buildFromPublisher(queryWorkflowOutputs(id), Option(id), List(workflowOutputsComponent))
  }

  def workflowLogs(id: WorkflowId)(implicit ec: ExecutionContext): IO[MetadataComponent] = {
    buildFromPublisher(queryLogs(id), Option(id))
  }
}
