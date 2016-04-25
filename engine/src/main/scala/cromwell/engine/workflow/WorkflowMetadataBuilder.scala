package cromwell.engine.workflow

import akka.actor.ActorSystem
import cromwell.core.WorkflowId
import cromwell.database.obj.WorkflowExecution
import cromwell.engine
import cromwell.engine._
import cromwell.engine.backend.{CallMetadata, WorkflowDescriptor}
import cromwell.engine.db._
import cromwell.engine.finalcall.FinalCall._
import cromwell.engine.workflow.WorkflowManagerActor.WorkflowNotFoundException
import cromwell.webservice._
import org.joda.time.DateTime
import spray.json._
import wdl4s._

import scala.concurrent.Future

object WorkflowMetadataBuilder {

  private type DBMap[+V] = Map[ExecutionDatabaseKey, V]

  // TODO: This assertion could be added to the db layer: "In the future I'll fail if the workflow doesn't exist"
  def assertWorkflowExistence(id: WorkflowId, workflowState: Option[WorkflowState]): Future[Unit] = {
    // Confirm the workflow exists by querying its state.  If no state is found the workflow doesn't exist.
    workflowState match {
      case None => Future.failed(new WorkflowNotFoundException(s"Workflow '$id' not found"))
      case _ => Future.successful(Unit)
    }
  }
  private def build(workflowDescriptor: WorkflowDescriptor,
                    execution: WorkflowExecution,
                    workflowOutputs: engine.WorkflowOutputs,
                    callMetadata: Map[FullyQualifiedName, Seq[CallMetadata]],
                    workflowFailures: Traversable[FailureEventEntry]):
  WorkflowMetadataResponse = {

    val startDate = new DateTime(execution.startDt)
    val endDate = execution.endDt map { new DateTime(_) }
    val workflowInputs = workflowDescriptor.sourceFiles.inputsJson.parseJson.asInstanceOf[JsObject]
    val failures = if (workflowFailures.isEmpty) None else Option(workflowFailures)

    WorkflowMetadataResponse(
      id = execution.workflowExecutionUuid.toString,
      workflowName = execution.name,
      status = execution.status,
      // We currently do not make a distinction between the submission and start dates of a workflow, but it's
      // possible at least theoretically that a workflow might not begin to execute immediately upon submission.
      submission = startDate,
      start = Option(startDate),
      end = endDate,
      inputs = workflowInputs,
      outputs = Option(workflowOutputs) map { _.mapToValues },
      calls = callMetadata,
      failures = failures.map(_.toSeq))
  }

  private def callFailuresMap(failureEvents: Seq[QualifiedFailureEventEntry]): DBMap[Seq[FailureEventEntry]] = {
    failureEvents filter { _.execution.isDefined } groupBy { _.execution } map {
      case (key, qualifiedEntries: Seq[QualifiedFailureEventEntry]) => key.get -> (qualifiedEntries map { _.dequalify })
    }
  }

  private def emptySeq[A] = Future.successful(Seq.empty[A])

  private def retrieveTraversable[A](runQuery: Boolean,
                                     queryResult: => Future[Traversable[A]],
                                     defaultResult: => Future[Traversable[A]] = emptySeq): Future[Traversable[A]] = {
    if (runQuery) queryResult else defaultResult
  }

  private def emptyMap[V] = Future.successful(Map.empty[ExecutionDatabaseKey, V])

  private def retrieveMap[V](runQuery: Boolean,
                             queryResult: => Future[DBMap[V]],
                             defaultResult: => Future[DBMap[V]] = emptyMap): Future[DBMap[V]] = {
    if (runQuery) queryResult else defaultResult
  }
}

class WorkflowMetadataBuilder(id: WorkflowId, parameters: WorkflowMetadataQueryParameters)(implicit actorSystem: ActorSystem) {

  private[this] implicit val ec = actorSystem.dispatcher

  import cromwell.engine.db.DataAccess.globalDataAccess

  import WorkflowMetadataBuilder._

  private def futureAssertWorkflowExistsByState = for {
    workflowState <- globalDataAccess.getWorkflowState(id)
    // TODO: This assertion could be added to the db layer: "In the future I'll fail if the workflow doesn't exist"
    _ <- assertWorkflowExistence(id, workflowState)
  } yield ()

  private def futureWorkflowExecutionAndAux = globalDataAccess.getWorkflowExecutionAndAux(id)

  // If outputs requested, don't retrieve the execution infos, only executions
  private def futureInfosByExecution = retrieveTraversable(parameters.outputs,
    globalDataAccess.infosByExecution(id),
    globalDataAccess.getExecutions(id) map { _ map { ExecutionInfosByExecution(_, Seq.empty) } })

  // If outputs requested, get the workflow outputs
  private def futureWorkflowOutputs = retrieveTraversable(parameters.outputs, globalDataAccess.getWorkflowOutputs(id))

  // If timings requested, we need the call statuses and execution events
  private def futureCallToStatusMap = retrieveMap(parameters.timings, globalDataAccess.getExecutionStatuses(id))

  private def futureExecutionEvents = retrieveMap(parameters.timings, globalDataAccess.getAllExecutionEvents(id))

  // Drop the below if timings _or_ outputs requested
  private def futureCallInputs = retrieveTraversable(parameters.timings && parameters.outputs,
    globalDataAccess.getAllInputs(id))

  private def futureCallOutputs = retrieveTraversable(parameters.timings && parameters.outputs,
    globalDataAccess.getAllOutputs(id))

  private def futureCallCacheData = retrieveTraversable(parameters.timings && parameters.outputs,
    globalDataAccess.callCacheDataByExecution(id))

  private def futureRuntimeAttributes = retrieveMap(parameters.timings && parameters.outputs,
    globalDataAccess.getAllRuntimeAttributes(id))

  private def futureFailures = retrieveTraversable(parameters.timings && parameters.outputs,
    globalDataAccess.getFailureEvents(id))

  def build(): Future[WorkflowMetadataResponse] = {

    for {
      assertWorkflowExistsByState <- futureAssertWorkflowExistsByState
      workflowExecutionAndAux <- futureWorkflowExecutionAndAux
      infosByExecution <- futureInfosByExecution
      workflowOutputs <- futureWorkflowOutputs
      callToStatusMap <- futureCallToStatusMap
      executionEvents <- futureExecutionEvents
      callInputs <- futureCallInputs
      callOutputs <- futureCallOutputs
      callCacheData <- futureCallCacheData
      runtimeAttributes <- futureRuntimeAttributes
      failures <- futureFailures

      // Database work complete, but we do need one more future to get the workflow descriptor
      workflowDescriptor <- workflowDescriptorFromExecutionAndAux(workflowExecutionAndAux)
    } yield {

      val execution = workflowExecutionAndAux.execution
      val nonFinalEvents = executionEvents.filterKeys(!_.fqn.isFinalCall)
      val nonFinalInfosByExecution = infosByExecution.filterNot(_.execution.callFqn.isFinalCall)

      val wfFailures = failures collect {
        case QualifiedFailureEventEntry(_, None, message, timestamp) => FailureEventEntry(message, timestamp)
      }
      val callFailures = callFailuresMap(failures.toSeq)

      val engineWorkflowOutputs = SymbolStoreEntry.toWorkflowOutputs(workflowOutputs)
      val callMetadata = CallMetadataBuilder.build(nonFinalInfosByExecution, callInputs, callOutputs, nonFinalEvents,
        runtimeAttributes, callCacheData, callFailures)

      WorkflowMetadataBuilder.build(workflowDescriptor, execution, engineWorkflowOutputs, callMetadata, wfFailures)
    }

  }
}
