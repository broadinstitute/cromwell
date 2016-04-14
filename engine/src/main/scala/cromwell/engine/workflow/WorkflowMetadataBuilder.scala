package cromwell.engine.workflow

import akka.actor.ActorSystem
import cromwell.core.WorkflowId
import cromwell.engine
import cromwell.engine._
import cromwell.engine.backend.CallMetadata
import cromwell.engine.db.DataAccess.WorkflowExecutionAndAux
import cromwell.engine.db.slick._
import cromwell.engine.db.{CallStatus, ExecutionDatabaseKey, ExecutionInfosByExecution}
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
                                    workflowExecution: WorkflowExecution,
                                    workflowOutputs: Traversable[SymbolStoreEntry],
                                    workflowExecutionAux: WorkflowExecutionAux,
                                    executionStatuses: Map[ExecutionDatabaseKey, CallStatus],
                                    callInputs: Traversable[SymbolStoreEntry],
                                    callOutputs: Traversable[SymbolStoreEntry],
                                    infosByExecution: Traversable[ExecutionInfosByExecution],
                                    runtimeAttributes: Map[ExecutionDatabaseKey, Map[String, String]],
                                    executionEvents: Map[ExecutionDatabaseKey, Seq[ExecutionEventEntry]],
                                    cacheData: Traversable[ExecutionWithCacheData],
                                    failures: Seq[QualifiedFailureEventEntry]) (implicit ec: ExecutionContext, actorSystem: ActorSystem):
  Future[WorkflowMetadataResponse] = {
    val executionAndAux = WorkflowExecutionAndAux(workflowExecution, workflowExecutionAux)
    workflowDescriptorFromExecutionAndAux(executionAndAux) map { workflowDescriptor =>
      val nonFinalRuntimeAttributes = runtimeAttributes.filterKeys(!_.isFinalCall)
      val nonFinalEvents = executionEvents.filterKeys(!_.isFinalCall)
      val nonFinalInfosByExecution = infosByExecution.filterNot(_.execution.toKey.isFinalCall)

      val wfFailures = failures collect {
        case QualifiedFailureEventEntry(_, None, message, timestamp) => FailureEventEntry(message, timestamp)
      }
      val callFailures = callFailuresMap(failures)

      val engineWorkflowOutputs = SymbolStoreEntry.toWorkflowOutputs(workflowOutputs)
      val callMetadata = CallMetadataBuilder.build(nonFinalInfosByExecution, callInputs, callOutputs, nonFinalEvents,
        nonFinalRuntimeAttributes, cacheData, callFailures)
      buildWorkflowMetadata(workflowExecution, workflowExecutionAux, engineWorkflowOutputs, callMetadata, wfFailures)
    }
  }

  def workflowMetadata(id: WorkflowId)(implicit ec: ExecutionContext, actorSystem: ActorSystem): Future[WorkflowMetadataResponse] = {

    // TODO: This entire block of chained database actions should be a single request to the db layer.
    import cromwell.engine.db.DataAccess.globalDataAccess

    for {
      workflowState <- globalDataAccess.getWorkflowState(id)
      // TODO: This assertion could be added to the db layer: "In the future I'll fail if the workflow doesn't exist"
      _ <- assertWorkflowExistence(id, workflowState)
      workflowExecution <- globalDataAccess.getWorkflowExecution(id)
      workflowOutputs <- globalDataAccess.getWorkflowOutputs(id)
      workflowExecutionAux <- globalDataAccess.getWorkflowExecutionAux(id)
      callToStatusMap <- globalDataAccess.getExecutionStatuses(id)
      callInputs <- globalDataAccess.getAllInputs(id)
      callOutputs <- globalDataAccess.getAllOutputs(id)
      infosByExecution <- globalDataAccess.infosByExecution(id)
      callCacheData <- globalDataAccess.callCacheDataByExecution(id)
      runtimeAttributes <- globalDataAccess.getAllRuntimeAttributes(id)
      executionEvents <- globalDataAccess.getAllExecutionEvents(id)
      failures <- globalDataAccess.getFailureEvents(id)
      wfMetadata <- buildWorkflowMetadata(id, workflowExecution, workflowOutputs, workflowExecutionAux,
        callToStatusMap, callInputs, callOutputs, infosByExecution, runtimeAttributes, executionEvents, callCacheData, failures)
    } yield wfMetadata
  }
}
