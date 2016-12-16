package cromwell.engine.workflow.lifecycle.execution

import java.time.OffsetDateTime

import akka.actor.ActorRef
import cromwell.backend.BackendJobDescriptorKey
import cromwell.core.ExecutionStatus._
import cromwell.core._
import cromwell.services.metadata.MetadataService._
import cromwell.services.metadata._
import wdl4s._
import wdl4s.values.WdlValue

import scala.util.Random

trait CallMetadataHelper {
  
  def workflowIdForCallMetadata: WorkflowId
  def serviceRegistryActor: ActorRef

  def pushNewCallMetadata(callKey: CallKey, backendName: Option[String]) = {
    val startEvents = List(
      Option(MetadataEvent(metadataKeyForCall(callKey, CallMetadataKeys.Start), MetadataValue(OffsetDateTime.now))),
      backendName map { name => MetadataEvent(metadataKeyForCall(callKey, CallMetadataKeys.Backend), MetadataValue(name)) }
    ).flatten

    serviceRegistryActor ! PutMetadataAction(startEvents)
  }

  def pushQueuedCallMetadata(diffs: Seq[WorkflowExecutionDiff]) = {
    val startingEvents = for {
      diff <- diffs
      (jobKey, executionState) <- diff.executionStoreChanges if jobKey.isInstanceOf[BackendJobDescriptorKey] && executionState == ExecutionStatus.QueuedInCromwell
    } yield MetadataEvent(metadataKeyForCall(jobKey, CallMetadataKeys.ExecutionStatus), MetadataValue(ExecutionStatus.QueuedInCromwell))
    serviceRegistryActor ! PutMetadataAction(startingEvents)
  }

  def pushStartingCallMetadata(callKey: CallKey) = {
    val statusChange = MetadataEvent(metadataKeyForCall(callKey, CallMetadataKeys.ExecutionStatus), MetadataValue(ExecutionStatus.Starting))
    serviceRegistryActor ! PutMetadataAction(statusChange)
  }

  def pushRunningCallMetadata(key: CallKey, evaluatedInputs: EvaluatedTaskInputs) = {
    val inputEvents = evaluatedInputs match {
      case empty if empty.isEmpty =>
        List(MetadataEvent.empty(metadataKeyForCall(key, s"${CallMetadataKeys.Inputs}")))
      case inputs =>
        inputs flatMap {
          case (inputName, inputValue) =>
            wdlValueToMetadataEvents(metadataKeyForCall(key, s"${CallMetadataKeys.Inputs}:${inputName.unqualifiedName}"), inputValue)
        }
    }

    val runningEvent = List(MetadataEvent(metadataKeyForCall(key, CallMetadataKeys.ExecutionStatus), MetadataValue(ExecutionStatus.Running)))

    serviceRegistryActor ! PutMetadataAction(runningEvent ++ inputEvents)
  }

  def pushWorkflowOutputMetadata(outputs: Map[LocallyQualifiedName, WdlValue]) = {
    val events = outputs match {
      case empty if empty.isEmpty => List(MetadataEvent.empty(MetadataKey(workflowIdForCallMetadata, None, WorkflowMetadataKeys.Outputs)))
      case _ => outputs flatMap {
        case (outputName, outputValue) =>
          wdlValueToMetadataEvents(MetadataKey(workflowIdForCallMetadata, None, s"${WorkflowMetadataKeys.Outputs}:$outputName"), outputValue)
      }
    }

    serviceRegistryActor ! PutMetadataAction(events)
  }

  def pushSuccessfulCallMetadata(jobKey: JobKey, returnCode: Option[Int], outputs: CallOutputs) = {
    val completionEvents = completedCallMetadataEvents(jobKey, ExecutionStatus.Done, returnCode)

    val outputEvents = outputs match {
      case empty if empty.isEmpty =>
        List(MetadataEvent.empty(metadataKeyForCall(jobKey, s"${CallMetadataKeys.Outputs}")))
      case _ =>
        outputs flatMap { case (lqn, outputValue) => wdlValueToMetadataEvents(metadataKeyForCall(jobKey, s"${CallMetadataKeys.Outputs}:$lqn"), outputValue.wdlValue) }
    }

    serviceRegistryActor ! PutMetadataAction(completionEvents ++ outputEvents)
  }

  def pushFailedCallMetadata(jobKey: JobKey, returnCode: Option[Int], failure: Throwable, retryableFailure: Boolean) = {
    val failedState = if (retryableFailure) ExecutionStatus.Preempted else ExecutionStatus.Failed
    val completionEvents = completedCallMetadataEvents(jobKey, failedState, returnCode)
    val retryableFailureEvent = MetadataEvent(metadataKeyForCall(jobKey, CallMetadataKeys.RetryableFailure), MetadataValue(retryableFailure))
    val failureEvents = throwableToMetadataEvents(metadataKeyForCall(jobKey, s"${CallMetadataKeys.Failures}[$randomNumberString]"), failure).+:(retryableFailureEvent)

    serviceRegistryActor ! PutMetadataAction(completionEvents ++ failureEvents)
  }

  def pushExecutionEventsToMetadataService(jobKey: JobKey, eventList: Seq[ExecutionEvent]) = {
    def metadataEvent(k: String, value: Any) = {
      val metadataValue = MetadataValue(value)
      val metadataKey = metadataKeyForCall(jobKey, k)
      MetadataEvent(metadataKey, metadataValue)
    }
    
    eventList.headOption foreach { firstEvent =>
      // The final event is only used as the book-end for the final pairing so the name is never actually used...
      val offset = firstEvent.offsetDateTime.getOffset
      val now = OffsetDateTime.now.withOffsetSameInstant(offset)
      val lastEvent = ExecutionEvent("!!Bring Back the Monarchy!!", now)
      val tailedEventList = eventList :+ lastEvent
      val events = tailedEventList.sliding(2) flatMap {
        case Seq(eventCurrent, eventNext) =>
          val eventKey = s"executionEvents[$randomNumberString]"
          List(
            metadataEvent(s"$eventKey:description", eventCurrent.name),
            metadataEvent(s"$eventKey:startTime", eventCurrent.offsetDateTime),
            metadataEvent(s"$eventKey:endTime", eventNext.offsetDateTime)
          )
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
  
  private def metadataKeyForCall(jobKey: JobKey, myKey: String) = MetadataKey(workflowIdForCallMetadata, Option(MetadataJobKey(jobKey.scope.fullyQualifiedName, jobKey.index, jobKey.attempt)), myKey)

  private def randomNumberString: String = Random.nextInt.toString.stripPrefix("-")
  
}
