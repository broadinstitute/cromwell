package cromwell.webservice.metadata

import cats.effect.IO
import com.typesafe.scalalogging.LazyLogging
import cromwell.core.{WorkflowId, WorkflowMetadataKeys}
import cromwell.services.MetadataServicesStore
import cromwell.services.metadata._
import cromwell.services.metadata.impl.MetadataDatabaseAccess
import cromwell.webservice.metadata.MetadataComponent._
import fs2.interop.reactivestreams._
import fs2.{Pipe, Stream}
import org.reactivestreams.Publisher

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

  val DefaultTimeout = 30.seconds

  def buildFromPublisher(eventsPublisher: Publisher[MetadataEvent], workflowId: Option[WorkflowId] = None, extraComponents: Seq[MetadataComponent] = Seq.empty)(implicit ec: ExecutionContext): IO[MetadataComponent] = {
    build(publisherToEventStream(eventsPublisher), toMetadataComponent, workflowId, extraComponents)
  }

  def build(eventStream: fs2.Stream[IO, MetadataEvent],
            eventToComponent: MetadataEvent => MetadataComponent,
            workflowId: Option[WorkflowId],
            extraComponents: Seq[MetadataComponent] = Seq.empty)(implicit ec: ExecutionContext): IO[MetadataComponent] = {
    implicit val timer = cats.effect.IO.timer(ec)
    implicit val contextShift = cats.effect.IO.contextShift(ec)
    implicit val concurrent = cats.effect.IO.ioConcurrentEffect(contextShift)

    eventStream
      .interruptWhen[IO]((Stream.sleep_[IO](DefaultTimeout) ++ Stream(true)) through logTimeout(workflowId))
      .map(eventToComponent)
      .append(fs2.Stream.emits(extraComponents))
      .compile
      .foldMonoid
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

  def workflowMetadataQuery(query: MetadataQuery, extraComponents: Seq[MetadataComponent] = Seq.empty)(implicit ec: ExecutionContext): IO[MetadataComponent] = {
    build(queryToEventStream(query), toMetadataComponent(query, ec), Option(query.workflowId), extraComponents)
  }

  def workflowOutputs(id: WorkflowId)(implicit ec: ExecutionContext): IO[MetadataComponent] = {
    buildFromPublisher(queryWorkflowOutputs(id), Option(id), List(workflowOutputsComponent))
  }

  def workflowLogs(id: WorkflowId)(implicit ec: ExecutionContext): IO[MetadataComponent] = {
    buildFromPublisher(queryLogs(id), Option(id))
  }
}
