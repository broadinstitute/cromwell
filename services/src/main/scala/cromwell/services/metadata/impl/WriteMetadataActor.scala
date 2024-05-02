package cromwell.services.metadata.impl

import akka.actor.{ActorLogging, ActorRef, Props}
import cats.data.NonEmptyVector
import cromwell.core.Dispatcher.ServiceDispatcher
import cromwell.core.Mailbox.PriorityMailbox
import cromwell.core.WorkflowId
import cromwell.core.instrumentation.InstrumentationPrefixes
import cromwell.services.metadata.{MetadataEvent, MetadataValue}
import cromwell.services.metadata.MetadataService._
import cromwell.services.metadata.impl.MetadataStatisticsRecorder.MetadataStatisticsRecorderSettings
import cromwell.services.{EnhancedBatchActor, MetadataServicesStore}
import wdl.util.StringUtil

import scala.concurrent.duration._
import scala.util.{Failure, Success}

class WriteMetadataActor(override val batchSize: Int,
                         override val flushRate: FiniteDuration,
                         override val serviceRegistryActor: ActorRef,
                         override val threshold: Int,
                         metadataStatisticsRecorderSettings: MetadataStatisticsRecorderSettings,
                         metadataKeysToClean: List[String]
) extends EnhancedBatchActor[MetadataWriteAction](flushRate, batchSize)
    with ActorLogging
    with MetadataDatabaseAccess
    with MetadataServicesStore {

  private val statsRecorder = MetadataStatisticsRecorder(metadataStatisticsRecorderSettings)

  override def process(e: NonEmptyVector[MetadataWriteAction]) = instrumentedProcess {
    val cleanedMetadataWriteActions = if (metadataKeysToClean.isEmpty) e else sanitizeInputs(e)
    val empty = (Vector.empty[MetadataEvent], List.empty[(Iterable[MetadataEvent], ActorRef)])

    val (putWithoutResponse, putWithResponse) = cleanedMetadataWriteActions.foldLeft(empty) {
      case ((putEvents, putAndRespondEvents), action: PutMetadataAction) =>
        (putEvents ++ action.events, putAndRespondEvents)
      case ((putEvents, putAndRespondEvents), action: PutMetadataActionAndRespond) =>
        (putEvents, putAndRespondEvents :+ (action.events -> action.replyTo))
    }
    val allPutEvents: Iterable[MetadataEvent] = putWithoutResponse ++ putWithResponse.flatMap(_._1)
    val dbAction = addMetadataEvents(allPutEvents)

    statsRecorder.processEventsAndGenerateAlerts(allPutEvents) foreach (a =>
      log.warning(s"${a.workflowId} has logged a heavy amount of metadata (${a.count} rows)")
    )

    dbAction onComplete {
      case Success(_) =>
        putWithResponse foreach { case (ev, replyTo) => replyTo ! MetadataWriteSuccess(ev) }
      case Failure(cause) =>
        val (outOfTries, stillGood) = cleanedMetadataWriteActions.toVector.partition(_.maxAttempts <= 1)

        handleOutOfTries(outOfTries, cause)
        handleEventsToReconsider(stillGood)
    }

    dbAction.map(_ => allPutEvents.size)
  }

  def sanitizeInputs(
    metadataWriteActions: NonEmptyVector[MetadataWriteAction]
  ): NonEmptyVector[MetadataWriteAction] =
    metadataWriteActions.map { metadataWriteAction =>
      val metadataEvents =
        metadataWriteAction.events.map { event =>
          event.value match {
            case Some(eventVal) => event.copy(value = Option(MetadataValue(StringUtil.cleanUtf8mb4(eventVal.value))))
            case None => event
          }
        }
      metadataWriteAction match {
        case action: PutMetadataAction => action.copy(events = metadataEvents)
        case actionAndResp: PutMetadataActionAndRespond => actionAndResp.copy(events = metadataEvents)
      }
    }

  private def countActionsByWorkflow(writeActions: Vector[MetadataWriteAction]): Map[WorkflowId, Int] =
    writeActions.flatMap(_.events).groupBy(_.key.workflowId).map { case (k, v) => k -> v.size }

  private def enumerateWorkflowWriteFailures(writeActions: Vector[MetadataWriteAction]): String =
    countActionsByWorkflow(writeActions).map { case (wfid, size) => s"$wfid: $size" }.mkString(", ")

  private def handleOutOfTries(writeActions: Vector[MetadataWriteAction], reason: Throwable): Unit = if (
    writeActions.nonEmpty
  ) {
    log.error(
      reason,
      "Metadata event writes have failed irretrievably for the following workflows. They will be lost: " + enumerateWorkflowWriteFailures(
        writeActions
      )
    )

    writeActions foreach {
      case PutMetadataActionAndRespond(ev, replyTo, _) => replyTo ! MetadataWriteFailure(reason, ev)
      case _: PutMetadataAction => () // We need to satisfy the exhaustive match but there's nothing special to do here
    }
  }

  private def handleEventsToReconsider(writeActions: Vector[MetadataWriteAction]): Unit = if (writeActions.nonEmpty) {
    log.warning(
      "Metadata event writes have failed for the following workflows. They will be retried: " + enumerateWorkflowWriteFailures(
        writeActions
      )
    )

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
  def commandToData(snd: ActorRef): PartialFunction[Any, MetadataWriteAction] = { case command: MetadataWriteAction =>
    command
  }
}

object WriteMetadataActor {

  def props(dbBatchSize: Int,
            flushRate: FiniteDuration,
            serviceRegistryActor: ActorRef,
            threshold: Int,
            statisticsRecorderSettings: MetadataStatisticsRecorderSettings,
            metadataKeysToClean: List[String]
  ): Props =
    Props(
      new WriteMetadataActor(dbBatchSize,
                             flushRate,
                             serviceRegistryActor,
                             threshold,
                             statisticsRecorderSettings,
                             metadataKeysToClean
      )
    )
      .withDispatcher(ServiceDispatcher)
      .withMailbox(PriorityMailbox)
}
