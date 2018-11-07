package cromwell.engine.workflow.lifecycle.execution

import java.time.OffsetDateTime

import akka.actor.ActorRef
import cromwell.backend.async.JobAlreadyFailedInJobStore
import cromwell.core.ExecutionStatus._
import cromwell.core._
import cromwell.services.metadata.MetadataService._
import cromwell.services.metadata._
import wdl.draft2.model._
import wom.values.{WomEvaluatedCallInputs, WomValue}

import scala.util.Random

trait CallMetadataHelper {
  
  def workflowIdForCallMetadata: WorkflowId
  def serviceRegistryActor: ActorRef

  def pushNewCallMetadata(callKey: CallKey, backendName: Option[String], serviceRegistryActor: ActorRef) = {
    val startEvents = List(
      Option(MetadataEvent(metadataKeyForCall(callKey, CallMetadataKeys.Start), MetadataValue(OffsetDateTime.now))),
      Option(MetadataEvent(metadataKeyForCall(callKey, CallMetadataKeys.ExecutionStatus), MetadataValue(ExecutionStatus.QueuedInCromwell))),
      backendName map { name => MetadataEvent(metadataKeyForCall(callKey, CallMetadataKeys.Backend), MetadataValue(name)) }
    ).flatten

    serviceRegistryActor ! PutMetadataAction(startEvents)
  }

  def pushStartingCallMetadata(callKey: CallKey) = {
    val statusChange = MetadataEvent(metadataKeyForCall(callKey, CallMetadataKeys.ExecutionStatus), MetadataValue(ExecutionStatus.Starting))
    serviceRegistryActor ! PutMetadataAction(statusChange)
  }

  def pushRunningCallMetadata(key: CallKey, evaluatedInputs: WomEvaluatedCallInputs) = {
    val inputEvents = evaluatedInputs match {
      case empty if empty.isEmpty =>
        List(MetadataEvent.empty(metadataKeyForCall(key, s"${CallMetadataKeys.Inputs}")))
      case inputs =>
        inputs flatMap {
          case (inputName, inputValue) =>
            womValueToMetadataEvents(metadataKeyForCall(key, s"${CallMetadataKeys.Inputs}:${inputName.name}"), inputValue)
        }
    }

    val runningEvent = List(MetadataEvent(metadataKeyForCall(key, CallMetadataKeys.ExecutionStatus), MetadataValue(ExecutionStatus.Running)))

    serviceRegistryActor ! PutMetadataAction(runningEvent ++ inputEvents)
  }

  def pushWorkflowOutputMetadata(outputs: Map[LocallyQualifiedName, WomValue]) = {
    val events = if (outputs.isEmpty) {
      List(MetadataEvent.empty(MetadataKey(workflowIdForCallMetadata, None, WorkflowMetadataKeys.Outputs)))
    } else {
      outputs flatMap { case (outputName, outputValue) =>
        womValueToMetadataEvents(MetadataKey(workflowIdForCallMetadata, None, s"${WorkflowMetadataKeys.Outputs}:$outputName"), outputValue)
      }
    }

    serviceRegistryActor ! PutMetadataAction(events)
  }
  
  def pushWaitingForQueueSpaceCallMetadata(jobKey: JobKey) = {
    val event = MetadataEvent(metadataKeyForCall(jobKey, CallMetadataKeys.ExecutionStatus), MetadataValue(WaitingForQueueSpace))
    serviceRegistryActor ! PutMetadataAction(event)
  }

  def pushSuccessfulCallMetadata(jobKey: JobKey, returnCode: Option[Int], outputs: CallOutputs) = {
    val completionEvents = completedCallMetadataEvents(jobKey, ExecutionStatus.Done, returnCode)

    val outputEvents = outputs.outputs match {
      case empty if empty.isEmpty =>
        List(MetadataEvent.empty(metadataKeyForCall(jobKey, s"${CallMetadataKeys.Outputs}")))
      case _ =>
        outputs.outputs flatMap { case (outputPort, outputValue) =>
          womValueToMetadataEvents(metadataKeyForCall(jobKey, s"${CallMetadataKeys.Outputs}:${outputPort.internalName}"), outputValue)
        }
    }

    serviceRegistryActor ! PutMetadataAction(completionEvents ++ outputEvents)
  }

  def pushFailedCallMetadata(jobKey: JobKey, returnCode: Option[Int], failure: Throwable, retryableFailure: Boolean) = {
    val failedState = if (retryableFailure) ExecutionStatus.RetryableFailure else ExecutionStatus.Failed
    val completionEvents = completedCallMetadataEvents(jobKey, failedState, returnCode)
    val retryableFailureEvent = MetadataEvent(metadataKeyForCall(jobKey, CallMetadataKeys.RetryableFailure), MetadataValue(retryableFailure))

    val failureEvents = failure match {
      // If the job was already failed, don't republish the failure reasons, they're already there
      case _: JobAlreadyFailedInJobStore => List.empty
      case _ =>
        throwableToMetadataEvents(metadataKeyForCall(jobKey, s"${CallMetadataKeys.Failures}"), failure).+:(retryableFailureEvent)
    }
  
    serviceRegistryActor ! PutMetadataAction(completionEvents ++ failureEvents)
  }

  def pushAbortedCallMetadata(jobKey: JobKey) = {
    val completionEvents = completedCallMetadataEvents(jobKey, ExecutionStatus.Aborted, None)

    serviceRegistryActor ! PutMetadataAction(completionEvents)
  }

  def pushExecutionEventsToMetadataService(jobKey: JobKey, eventList: Seq[ExecutionEvent]) = {
    def metadataEvent(k: String, value: Any) = {
      val metadataValue = MetadataValue(value)
      val metadataKey = metadataKeyForCall(jobKey, k)
      MetadataEvent(metadataKey, metadataValue)
    }
    
    val sortedEvents = eventList.sortBy(_.offsetDateTime)

    sortedEvents.headOption foreach { firstEvent =>
      // The final event is only used as the book-end for the final pairing so the name is never actually used...
      val offset = firstEvent.offsetDateTime.getOffset
      val now = OffsetDateTime.now.withOffsetSameInstant(offset)
      val lastEvent = ExecutionEvent("!!Bring Back the Monarchy!!", now)
      val tailedEventList = sortedEvents :+ lastEvent
      val events = tailedEventList.sliding(2) flatMap {
        case Seq(eventCurrent, eventNext) =>
          val eventKey = s"${CallMetadataKeys.ExecutionEvents}[$randomNumberString]"
          List(
            metadataEvent(s"$eventKey:description", eventCurrent.name),
            metadataEvent(s"$eventKey:startTime", eventCurrent.offsetDateTime),
            metadataEvent(s"$eventKey:endTime", eventNext.offsetDateTime)
          ) ++ (eventCurrent.grouping map { g => metadataEvent(s"$eventKey:grouping", g) })
      }

      serviceRegistryActor ! PutMetadataAction(events.toIterable)
    }
  }

  private def completedCallMetadataEvents(jobKey: JobKey, executionStatus: ExecutionStatus, returnCode: Option[Int]) = {
    val returnCodeEvent = returnCode map { rc =>
      List(MetadataEvent(metadataKeyForCall(jobKey, CallMetadataKeys.ReturnCode), MetadataValue(rc)))
    }

    List(
      MetadataEvent(metadataKeyForCall(jobKey, CallMetadataKeys.ExecutionStatus), MetadataValue(executionStatus)),
      MetadataEvent(metadataKeyForCall(jobKey, CallMetadataKeys.End), MetadataValue(OffsetDateTime.now))
    ) ++ returnCodeEvent.getOrElse(List.empty)
  }
  
  private def randomNumberString: String = Random.nextInt.toString.stripPrefix("-")

  def metadataKeyForCall(jobKey: JobKey, myKey: String) = MetadataKey(workflowIdForCallMetadata, Option(MetadataJobKey(jobKey.node.fullyQualifiedName, jobKey.index, jobKey.attempt)), myKey)
}
