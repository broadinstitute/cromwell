package cromwell.engine.workflow

import cromwell.engine
import cromwell.engine.ExecutionIndex._
import cromwell.engine._
import cromwell.engine.backend.{BackendCallJobDescriptor, Backend, CallMetadata, WorkflowLogs}
import cromwell.engine.db.slick._
import cromwell.engine.db.{CallStatus, ExecutionDatabaseKey}
import cromwell.engine.finalcall.FinalCall._
import cromwell.engine.workflow.WorkflowManagerActor._
import cromwell.webservice._
import org.joda.time.DateTime
import spray.json._
import wdl4s._

import scala.concurrent.{ExecutionContext, Future}
import scala.io.Source

object WorkflowMetadataBuilder {

  // TODO: This assertion could be added to the db layer: "In the future I'll fail if the workflow doesn't exist"
  def assertWorkflowExistence(id: WorkflowId, workflowState: Option[WorkflowState]): Future[Unit] = {
    // Confirm the workflow exists by querying its state.  If no state is found the workflow doesn't exist.
    workflowState match {
      case None => Future.failed(new WorkflowNotFoundException(s"Workflow '$id' not found"))
      case _ => Future.successful(Unit)
    }
  }

  private def assertCallFqnWellFormed(descriptor: WorkflowDescriptor, callFqn: FullyQualifiedName): String = {
    descriptor.namespace.resolve(callFqn) match {
      case Some(c: Call) => c.unqualifiedName
      case _ => throw new UnsupportedOperationException(
        s"Expected a fully qualified name to have at least two parts but got $callFqn")
    }
  }

  private def hasLogs(entries: Iterable[ExecutionDatabaseKey])(key: ExecutionDatabaseKey) = {
    !key.fqn.isScatter && !key.isCollector(entries)
  }

  def workflowStdoutStderr(workflowId: WorkflowId, backend: Backend,
                           workflowDescriptor: WorkflowDescriptor,
                           executionStatuses: Map[ExecutionDatabaseKey, CallStatus]): WorkflowLogs = {
    val statusMap = executionStatuses filterKeys {
      !_.fqn.isFinalCall
    } mapValues {
      _.executionStatus
    }

    val sortedMap = statusMap.toSeq.sortBy(_._1.index)
    val callsToPaths = for {
      (key, status) <- sortedMap if !key.fqn.isFinalCall && hasLogs(statusMap.keys)(key)
      callName = assertCallFqnWellFormed(workflowDescriptor, key.fqn)
      call = workflowDescriptor.namespace.workflow.findCallByName(callName).get
      callKey = BackendCallKey(call, key.index, key.attempt)
      callStandardOutput = backend.stdoutStderr(BackendCallJobDescriptor(workflowDescriptor, callKey))
    } yield key -> callStandardOutput

    /* Some FP "magic" to transform the pairs of (key, logs) into the final result:
    grouped by FQNS, ordered by shards, and then ordered by attempts */
    callsToPaths groupBy {
      _._1.fqn
    } mapValues {
      key => key.groupBy(_._1.index).values.toIndexedSeq.sortBy(_.head._1.index) map {
        _.sortBy(_._1.attempt).map(_._2).toIndexedSeq
      }
    }
  }

  private def buildWorkflowMetadata(workflowExecution: WorkflowExecution,
                                    workflowExecutionAux: WorkflowExecutionAux,
                                    workflowOutputs: engine.WorkflowOutputs,
                                    callMetadata: Map[FullyQualifiedName, Seq[CallMetadata]],
                                    workflowFailures: Seq[FailureEventEntry]):
  WorkflowMetadataResponse = {

    val startDate = new DateTime(workflowExecution.startDt)
    val endDate = workflowExecution.endDt map {
      new DateTime(_)
    }
    // TODO: Refactor db to expose the clob <-> string implicits.
    val workflowInputs = Source.fromInputStream(workflowExecutionAux.jsonInputs.getAsciiStream)
      .mkString.parseJson.asInstanceOf[JsObject]
    val failures = if (workflowFailures.isEmpty) None else Option(workflowFailures)

    WorkflowMetadataResponse(
      id = workflowExecution.workflowExecutionUuid.toString,
      workflowName = workflowExecution.name,
      status = workflowExecution.status,
      // We currently do not make a distinction between the submission and start dates of a workflow, but it's
      // possible at least theoretically that a workflow might not begin to execute immediately upon submission.
      submission = startDate,
      start = Option(startDate),
      end = endDate,
      inputs = workflowInputs,
      outputs = Option(workflowOutputs) map {
        _.mapToValues
      },
      calls = callMetadata,
      failures = failures)
  }

  private def callFailuresMap(failureEvents: Seq[QualifiedFailureEventEntry]): Map[ExecutionDatabaseKey, Seq[FailureEventEntry]] = {
    failureEvents filter { _.execution.isDefined } groupBy { _.execution } map {
      case (key, qualifiedEntries: Seq[QualifiedFailureEventEntry]) => key.get -> (qualifiedEntries map { _.dequalify })
    }
  }

  private def buildWorkflowMetadata(id: WorkflowId,
                                    backend: Backend,
                                    workflowExecution: WorkflowExecution,
                                    workflowOutputs: Traversable[SymbolStoreEntry],
                                    workflowExecutionAux: WorkflowExecutionAux,
                                    workflowDescriptor: WorkflowDescriptor,
                                    executionStatuses: Map[ExecutionDatabaseKey, CallStatus],
                                    callInputs: Traversable[SymbolStoreEntry],
                                    callOutputs: Traversable[SymbolStoreEntry],
                                    infosByExecution: Traversable[ExecutionInfosByExecution],
                                    runtimeAttributes: Map[ExecutionDatabaseKey, Map[String, String]],
                                    executionEvents: Map[ExecutionDatabaseKey, Seq[ExecutionEventEntry]],
                                    failures: Seq[QualifiedFailureEventEntry]):
  WorkflowMetadataResponse = {
    val nonFinalEvents = executionEvents.filterKeys(!_.fqn.isFinalCall)
    val nonFinalStatuses = executionStatuses.filterKeys(!_.fqn.isFinalCall)
    val nonFinalInfosByExecution = infosByExecution.filterNot(_.execution.callFqn.isFinalCall)

    val wfFailures = failures collect {
      case QualifiedFailureEventEntry(_, None, message, timestamp) => FailureEventEntry(message, timestamp)
    }
    val callFailures = callFailuresMap(failures)

    val engineWorkflowOutputs = SymbolStoreEntry.toWorkflowOutputs(workflowOutputs)
    val callStandardStreamsMap = workflowStdoutStderr(id, backend, workflowDescriptor, nonFinalStatuses)
    val callMetadata = CallMetadataBuilder.build(nonFinalInfosByExecution, callStandardStreamsMap, callInputs,
      callOutputs, nonFinalEvents, runtimeAttributes, callFailures)
    buildWorkflowMetadata(workflowExecution, workflowExecutionAux, engineWorkflowOutputs, callMetadata, wfFailures)
  }

  def workflowMetadata(id: WorkflowId, backend: Backend)
                      (implicit ec: ExecutionContext): Future[WorkflowMetadataResponse] = {

    // TODO: This entire block of chained database actions should be a single request to the db layer.
    import cromwell.engine.db.DataAccess.globalDataAccess

    for {
      workflowState <- globalDataAccess.getWorkflowState(id)
      // TODO: This assertion could be added to the db layer: "In the future I'll fail if the workflow doesn't exist"
      _ <- assertWorkflowExistence(id, workflowState)
      workflowExecution <- globalDataAccess.getWorkflowExecution(id)
      workflowOutputs <- globalDataAccess.getWorkflowOutputs(id)
      workflowExecutionAux <- globalDataAccess.getWorkflowExecutionAux(id)
      workflowDescriptor <- globalDataAccess.getWorkflow(id)
      callToStatusMap <- globalDataAccess.getExecutionStatuses(id)
      callInputs <- globalDataAccess.getAllInputs(id)
      callOutputs <- globalDataAccess.getAllOutputs(id)
      infosByExecution <- globalDataAccess.infosByExecution(id)
      runtimeAttributes <- globalDataAccess.getAllRuntimeAttributes(id)
      executionEvents <- globalDataAccess.getAllExecutionEvents(id)
      failures <- globalDataAccess.getFailureEvents(id)
    } yield buildWorkflowMetadata(
      id,
      backend,
      workflowExecution,
      workflowOutputs,
      workflowExecutionAux,
      workflowDescriptor,
      callToStatusMap,
      callInputs,
      callOutputs,
      infosByExecution,
      runtimeAttributes,
      executionEvents,
      failures)
  }
}
