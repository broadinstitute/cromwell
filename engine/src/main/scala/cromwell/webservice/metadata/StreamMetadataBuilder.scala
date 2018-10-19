package cromwell.webservice.metadata

import cats.effect.IO
import cats.effect.IO._
import cats.instances.tuple._
import cats.instances.vector._
import cats.syntax.parallel._
import com.typesafe.scalalogging.LazyLogging
import common.validation.Validation._
import cromwell.core.WorkflowId
import cromwell.services.MetadataServicesStore
import cromwell.services.metadata.MetadataService.GetSingleWorkflowMetadataAction
import cromwell.services.metadata._
import cromwell.services.metadata.impl.StreamMetadataDatabaseAccess
import cromwell.webservice.metadata.MetadataComponent._
import fs2.interop.reactivestreams._
import fs2.{Pipe, Stream => FStream}
import org.reactivestreams.Publisher
import spray.json._

import scala.concurrent.ExecutionContext
import scala.concurrent.duration._

/**
  * Builds a metadata JSON, streaming elements from the Database.
  * DO NOT stream the json itself back to the caller.
  */
class StreamMetadataBuilder(requestTimeout: FiniteDuration) extends StreamMetadataDatabaseAccess with MetadataServicesStore with LazyLogging {
  // A MetadataComponent for the workflow id
  private def workflowIdComponent(id: WorkflowId): MetadataComponent = MetadataObject("id" -> MetadataPrimitive(MetadataValue(id.toString)))
  
  // An empty calls section
  private val emptyCallSection: MetadataComponent = MetadataObject(CallsKey -> MetadataEmptyComponent)

  /**
    * Builds the json corresponding to the provided MetadataQuery
    */
  def workflowMetadataQuery(query: MetadataQuery)(implicit ec: ExecutionContext): IO[JsObject] = {
    buildQuery(query, List(workflowIdComponent(query.workflowId), emptyCallSection))
  }

  /**
    * Builds the json corresponding to the provided GetSingleWorkflowMetadataAction
    */
  def workflowMetadataQuery(action: GetSingleWorkflowMetadataAction)(implicit ec: ExecutionContext): IO[JsObject] = action match {
    case GetSingleWorkflowMetadataAction(workflowId, includeKeysOption, excludeKeysOption, expandSubWorkflows) =>
      val includeKeys = if (expandSubWorkflows) includeKeysOption map { _.::(CallMetadataKeys.SubWorkflowId) } else includeKeysOption
      workflowMetadataQuery(MetadataQuery(workflowId, None, None, includeKeys, excludeKeysOption, expandSubWorkflows))
  }

  /**
    * Builds the JsObject corresponding to a MetadataQuery, in a streaming fashion.
    * Uses a reactive stream publisher that gets turned into an fs2 stream which we use to pull the metadata rows from the database
    * and build a MetadataComponent as the rows get pulled. Use the MetadataComponent Monoid to fold them into one.
    * Ultimately serializes it to Json.
    *
    * Note: Unfortunately we cannot simply flatMap the sub workflow id streams into the main one.
    * That is because the database connection providing the stream will stay open and "locked" until the stream is consumed entirely (or cancelled).
    * Because of that if we make streams depend on each other we can create deadlocks.
    * Instead, this accumulates the sub workflow IDs we encounter in a Vector. Once the main stream is complete, we build the json corresponding
    * to those subworkflows recursively and then use them when serializing the main metadata from MetadataComponent to Json.
    */
  private def buildQuery(query: MetadataQuery, extraComponents: Seq[MetadataComponent])
                        (implicit ec: ExecutionContext): IO[JsObject] = {
    implicit val contextShift = cats.effect.IO.contextShift(ec)

    // Given a vector of workflow ids, build the corresponding jsons in parallel using the same query
    def buildSubWorkflows(ids: Vector[WorkflowId]): IO[Map[WorkflowId, JsObject]] = {
      ids.parTraverse[IO, IO.Par, (WorkflowId, JsObject)](id =>
        buildQuery(query.copy(workflowId = id), extraComponents).map(id -> _)
      ).map(_.toMap)
    }

    // Transform a MetadataEvent to a MetadataComponent, and optionally add it to the sub workflow id list if necessary
    def mapEvent(subWorkflowIds: Vector[WorkflowId], event: MetadataEvent): (Vector[WorkflowId], MetadataComponent) = {
      if (event.isSubWorkflowId)
        (event.extractSubWorkflowId.map(subWorkflowIds.:+(_)).getOrElse(subWorkflowIds), toMetadataComponent(event))
      else{
        val a = (subWorkflowIds, toMetadataComponent(event))
        a
      }
    }

    // Run the stream by mapping metadata events to metadata components while accumulating sub workflow ids. Then fold everything together
    def runStream(stream: FStream[IO, MetadataEvent]): IO[(Vector[WorkflowId], MetadataComponent)] = stream
      .mapAccumulate[Vector[WorkflowId], MetadataComponent](Vector.empty[WorkflowId])(mapEvent)
      .compile
      .foldMonoid
    
    def addExtraComponents(component: MetadataComponent) = component match {
      // If the stream was empty, result should be empty json, no extras
      case MetadataEmptyComponent => component
      case _ => extraComponents.toVector.fold(component)(MetadataComponentMonoid.combine)
    }

    for {
      // Get the publisher of metadata events from the database
      publisher <- queryToPublisher(query)
      // Turn it into a stream
      stream <- publisherToEventStream(publisher, Option(query.workflowId))
      // Read the stream and turn all that into a MetadataComponent along with the sub workflow ids in a vector
      subWorkflowIdsAndComponent <- runStream(stream)
      // Just unapply the Tuple for clarity
      (subWorkflowIds, component) = subWorkflowIdsAndComponent
      // Add extra components if any
      withExtras = addExtraComponents(component)
      // This is the Vector[WorkflowId] we recursively turn into its own Map[WorkflowId, JsObject]
      subWorkflowJsons <- buildSubWorkflows(subWorkflowIds)
      // Now we have everything, create a MetadataComponent Json serializer using the subworkflow jsons
      componentWriter = new MetadataComponentJsonWriter(subWorkflowJsons)
    } yield withExtras.toJson(componentWriter).asJsObject
  }

  /**
    * Obtain the Publisher from the database for a given query
    */
  private [metadata] def queryToPublisher(query: MetadataQuery): IO[Publisher[MetadataEvent]] = {
    queryMetadataEvents(query).toIO("Invalid metadata query")
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
