package cromwell.services.metadata.impl

import akka.actor.{ActorLogging, ActorRef, Props}
import cats.data.NonEmptyVector
import com.google.common.cache.CacheBuilder
import cromwell.core.Dispatcher.ServiceDispatcher
import cromwell.core.Mailbox.PriorityMailbox
import cromwell.core.WorkflowId
import cromwell.core.instrumentation.InstrumentationPrefixes
import cromwell.services.metadata.{MetadataEvent, MetadataKey, MetadataString, MetadataValue}
import cromwell.services.metadata.MetadataService._
import cromwell.services.{EnhancedBatchActor, MetadataServicesStore}

import scala.concurrent.duration._
import java.time.{Duration => JDuration}
import java.util.UUID
import java.util.concurrent.Callable

import cromwell.services.metadata.impl.WriteMetadataActor.MetadataStatisticsRecorder

import scala.util.{Failure, Success, Try}


class WriteMetadataActor(override val batchSize: Int,
                         override val flushRate: FiniteDuration,
                         override val serviceRegistryActor: ActorRef,
                         override val threshold: Int)
  extends EnhancedBatchActor[MetadataWriteAction](flushRate, batchSize)
    with ActorLogging
    with MetadataDatabaseAccess
    with MetadataServicesStore {

  private val statsRecorder = new MetadataStatisticsRecorder()

  override def process(e: NonEmptyVector[MetadataWriteAction]) = instrumentedProcess {
    val empty = (Vector.empty[MetadataEvent], List.empty[(Iterable[MetadataEvent], ActorRef)])

    val (putWithoutResponse, putWithResponse) = e.foldLeft(empty)({
      case ((putEvents, putAndRespondEvents), action: PutMetadataAction) =>
        (putEvents ++ action.events, putAndRespondEvents)
      case ((putEvents, putAndRespondEvents), action: PutMetadataActionAndRespond) =>
        (putEvents, putAndRespondEvents :+ (action.events -> action.replyTo))
    })
    val allPutEvents: Iterable[MetadataEvent] = putWithoutResponse ++ putWithResponse.flatMap(_._1)
    val dbAction = addMetadataEvents(allPutEvents)
    allPutEvents.groupBy(_.key.workflowId).foreach { case (id, list) =>
      val alerts = statsRecorder.recordNewRows(id, list)
      alerts.foreach(a => log.warning(s"${a.workflowId} has logged a heavy amount of metadata (${a.count} rows)"))
    }

    dbAction onComplete {
      case Success(_) =>
        putWithResponse foreach { case (ev, replyTo) => replyTo ! MetadataWriteSuccess(ev) }
      case Failure(cause) =>

        val (outOfTries, stillGood) = e.toVector.partition(_.maxAttempts <= 1)

        handleOutOfTries(outOfTries, cause)
        handleEventsToReconsider(stillGood)
    }

    dbAction.map(_ => allPutEvents.size)
  }

  private def countActionsByWorkflow(writeActions: Vector[MetadataWriteAction]): Map[WorkflowId, Int] =
    writeActions.flatMap(_.events).groupBy(_.key.workflowId).map { case (k, v) => k -> v.size }

  private def enumerateWorkflowWriteFailures(writeActions: Vector[MetadataWriteAction]): String =
    countActionsByWorkflow(writeActions).map { case (wfid, size) => s"$wfid: $size" }.mkString(", ")

  private def handleOutOfTries(writeActions: Vector[MetadataWriteAction], reason: Throwable): Unit = if (writeActions.nonEmpty) {
    log.error(reason, "Metadata event writes have failed irretrievably for the following workflows. They will be lost: " + enumerateWorkflowWriteFailures(writeActions))

    writeActions foreach {
      case PutMetadataActionAndRespond(ev, replyTo, _) => replyTo ! MetadataWriteFailure(reason, ev)
      case _: PutMetadataAction => () // We need to satisfy the exhaustive match but there's nothing special to do here
    }
  }

  private def handleEventsToReconsider(writeActions: Vector[MetadataWriteAction]): Unit = if (writeActions.nonEmpty) {
    log.warning("Metadata event writes have failed for the following workflows. They will be retried: " + enumerateWorkflowWriteFailures(writeActions))

    writeActions foreach {
      case action: PutMetadataAction => self ! action.copy(maxAttempts = action.maxAttempts - 1)
      case action: PutMetadataActionAndRespond => self ! action.copy(maxAttempts = action.maxAttempts - 1)
    }
  }

  // EnhancedBatchActor overrides
  override def receive = enhancedReceive.orElse(super.receive)
  override protected def weightFunction(command: MetadataWriteAction) = command.size
  override protected def instrumentationPath = MetadataServiceActor.MetadataInstrumentationPrefix
  override protected def instrumentationPrefix = InstrumentationPrefixes.ServicesPrefix
  def commandToData(snd: ActorRef): PartialFunction[Any, MetadataWriteAction] = {
    case command: MetadataWriteAction => command
  }
}

object WriteMetadataActor {

  def props(dbBatchSize: Int,
            flushRate: FiniteDuration,
            serviceRegistryActor: ActorRef,
            threshold: Int): Props =
    Props(new WriteMetadataActor(dbBatchSize, flushRate, serviceRegistryActor, threshold))
      .withDispatcher(ServiceDispatcher)
      .withMailbox(PriorityMailbox)

  final case class WorkflowMetadataWriteStatistics(totalWrites: Long, lastLogged: Long, knownParent: Option[WorkflowId])
  final case class HeavyMetadataAlert(workflowId: WorkflowId, count: Long)

  final class MetadataStatisticsRecorder() {

    val workflowCacheSize = 10L * 1000000 // TODO: Replace with the configured max workflows value (or max workflows x 2?)
    val metadataAlertInterval = 10L // TODO: Replace with our current understanding of a good metadata size

    // Statistics for each workflow
    private val metadataWriteStatisticsCache = CacheBuilder.newBuilder()
      .concurrencyLevel(2)
      .expireAfterAccess(JDuration.ofSeconds(4.hours.toSeconds))
      .maximumSize(workflowCacheSize)
      .build[WorkflowId, WorkflowMetadataWriteStatistics]()

    val writeStatisticsLoader: Callable[WorkflowMetadataWriteStatistics] = () => WorkflowMetadataWriteStatistics(0L, 0L, None)

    def recordNewRows(workflowId: WorkflowId, events: Iterable[MetadataEvent], fromSubworkflow: Boolean = false): Vector[HeavyMetadataAlert] = {
      val workflowWriteStats = metadataWriteStatisticsCache.get(workflowId, writeStatisticsLoader)
      val writesForWorkflow = workflowWriteStats.totalWrites + events.size.longValue()
      val knownParent: Option[WorkflowId] = workflowWriteStats.knownParent.orElse( if(fromSubworkflow) None else events.collectFirst {
        case MetadataEvent(MetadataKey(_, None, "parentWorkflowId"), Some(MetadataValue(value, MetadataString)), _) => Try(WorkflowId(UUID.fromString(value))).toOption
      }.flatten )

      knownParent.foreach(p => recordNewRows(p, events))

      if (writesForWorkflow > workflowWriteStats.lastLogged + metadataAlertInterval) {
        metadataWriteStatisticsCache.put(workflowId, workflowWriteStats.copy(totalWrites = writesForWorkflow, lastLogged = writesForWorkflow, knownParent = knownParent))
        Vector(HeavyMetadataAlert(workflowId, writesForWorkflow))
      } else {
        metadataWriteStatisticsCache.put(workflowId, workflowWriteStats.copy(totalWrites = writesForWorkflow, knownParent = knownParent))
        Vector.empty
      }
    }

    def recordSubworkflowRows(workflowId: WorkflowId, count: Long) = {

    }
  }
}
