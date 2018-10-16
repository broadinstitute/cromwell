package cromwell.webservice.metadata

import cats.effect.IO
import com.typesafe.scalalogging.StrictLogging
import cromwell.core.{WorkflowId, WorkflowMetadataKeys}
import cromwell.services.MetadataServicesStore
import cromwell.services.metadata.impl.MetadataDatabaseAccess
import cromwell.services.metadata.{MetadataEvent, MetadataValue, _}
import cromwell.webservice.metadata.MetadataComponent._
import fs2.interop.reactivestreams._
import org.reactivestreams.Publisher

import scala.concurrent.ExecutionContext
import scala.concurrent.duration._
object StreamMetadataBuilder extends StrictLogging with MetadataDatabaseAccess with MetadataServicesStore {
  
  val DefaultTimeout = 30.seconds
  
  def buildFromPublisher(eventsPublisher: Publisher[MetadataEvent], extraComponents: Seq[MetadataComponent] = Seq.empty)(implicit ec: ExecutionContext): IO[MetadataComponent] = {
    build(publisherToEventStream(eventsPublisher), toMetadataComponent, extraComponents)
  }

  def build(eventStream: fs2.Stream[IO, MetadataEvent],
            eventToComponent: MetadataEvent => MetadataComponent,
            extraComponents: Seq[MetadataComponent] = Seq.empty)(implicit ec: ExecutionContext): IO[MetadataComponent] = {
    implicit val timer = cats.effect.IO.timer(ec)
    val contextShift = cats.effect.IO.contextShift(ec)
    implicit val concurrent = cats.effect.IO.ioConcurrentEffect(contextShift)

    eventStream
      .interruptAfter(DefaultTimeout)
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
    build(queryToEventStream(query), toMetadataComponent(query, ec), extraComponents)
  }
  
  def workflowOutputs(id: WorkflowId)(implicit ec: ExecutionContext): IO[MetadataComponent] = {
    buildFromPublisher(queryWorkflowOutputs(id), List(workflowOutputsComponent))
  }

  def workflowLogs(id: WorkflowId)(implicit ec: ExecutionContext): IO[MetadataComponent] = {
    buildFromPublisher(queryLogs(id))
  }
}
