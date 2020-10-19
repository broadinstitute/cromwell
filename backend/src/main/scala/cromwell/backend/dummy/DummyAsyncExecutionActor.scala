package cromwell.backend.dummy

import java.time.OffsetDateTime
import java.util.UUID

import cats.data.Validated.{Invalid, Valid}
import cats.implicits._
import common.exception.AggregatedMessageException
import common.validation.ErrorOr.ErrorOr
import cromwell.backend.BackendJobLifecycleActor
import cromwell.backend.async.{ExecutionHandle, FailedNonRetryableExecutionHandle, PendingExecutionHandle, SuccessfulExecutionHandle}
import cromwell.backend.standard.{StandardAsyncExecutionActor, StandardAsyncExecutionActorParams, StandardAsyncJob}
import cromwell.core.CallOutputs
import cromwell.core.retry.SimpleExponentialBackoff
import wom.expression.NoIoFunctionSet
import wom.graph.GraphNodePort.{ExpressionBasedOutputPort, OutputPort}
import wom.values.WomValue

import scala.concurrent.Future
import scala.concurrent.duration._

class DummyAsyncExecutionActor(override val standardParams: StandardAsyncExecutionActorParams)
  extends BackendJobLifecycleActor
  with StandardAsyncExecutionActor  {

  /** The type of the run info when a job is started. */
  override type StandardAsyncRunInfo = String
  /** The type of the run status returned during each poll. */
  override type StandardAsyncRunState = String

  /** Should return true if the status contained in `thiz` is equivalent to `that`, delta any other data that might be carried around
    * in the state type.
    */
  override def statusEquivalentTo(thiz: String)(that: String): Boolean = thiz == that

  /**
    * Returns true when a job is complete, either successfully or unsuccessfully.
    *
    * @param runStatus The run status.
    * @return True if the job has completed.
    */
  override def isTerminal(runStatus: String): Boolean = runStatus == "DummyDone"

  override def dockerImageUsed: Option[String] = None

  override def pollBackOff: SimpleExponentialBackoff = SimpleExponentialBackoff(initialInterval = 1.second, maxInterval = 300.seconds, multiplier = 1.1)

  override def executeOrRecoverBackOff: SimpleExponentialBackoff = SimpleExponentialBackoff(initialInterval = 1.second, maxInterval = 300.seconds, multiplier = 1.1)


  var finishTime: Option[OffsetDateTime] = None

  override def executeAsync(): Future[ExecutionHandle] = {
    finishTime = Option(OffsetDateTime.now().plusMinutes(10))
    Future.successful(
      PendingExecutionHandle[StandardAsyncJob, StandardAsyncRunInfo, StandardAsyncRunState](
        jobDescriptor = jobDescriptor,
        pendingJob = StandardAsyncJob(UUID.randomUUID().toString),
        runInfo = Option("pending"),
        previousState = None
      )
    )
  }

  override def pollStatusAsync(handle: StandardAsyncPendingExecutionHandle): Future[String] = {
    finishTime match {
      case Some(ft) if (ft.isBefore(OffsetDateTime.now)) => Future.successful("done")
      case Some(_) => Future.successful("running")
      case None => Future.failed(new Exception("Dummy backend polling for status before finishTime is established(!!?)"))
    }

  }

  override def handlePollSuccess(oldHandle: StandardAsyncPendingExecutionHandle, state: String): Future[ExecutionHandle] = {

    if (state == "done") {
      val outputsValidation: ErrorOr[Map[OutputPort, WomValue]] = jobDescriptor.taskCall.outputPorts.toList.traverse {
        case expressionBasedOutputPort: ExpressionBasedOutputPort =>
          expressionBasedOutputPort.expression.evaluateValue(Map.empty, NoIoFunctionSet).map(expressionBasedOutputPort -> _)
        case other => s"Unknown output port type for Dummy backend output evaluator: ${other.getClass.getSimpleName}".invalidNel
      }.map(_.toMap)

      outputsValidation match {
        case Valid(outputs) =>
          Future.successful(SuccessfulExecutionHandle(
            outputs = CallOutputs(outputs.toMap),
            returnCode = 0,
            jobDetritusFiles = Map.empty,
            executionEvents = Seq.empty,
            resultsClonedFrom = None
          ))
        case Invalid(errors) =>
          Future.successful(FailedNonRetryableExecutionHandle(
            throwable = AggregatedMessageException("Evaluate outputs from dummy job", errors.toList),
            returnCode = None,
            kvPairsToSave = None
          ))
      }
    }
    else if (state == "running") {
      Future.successful(
        PendingExecutionHandle[StandardAsyncJob, StandardAsyncRunInfo, StandardAsyncRunState](
          jobDescriptor = jobDescriptor,
          pendingJob = StandardAsyncJob(UUID.randomUUID().toString),
          runInfo = Option("pending"),
          previousState = Option(state)
        )
      )
    }
    else {
      Future.failed(new Exception(s"Unexpected Dummy state in handlePollSuccess: $state"))
    }
  }
}
