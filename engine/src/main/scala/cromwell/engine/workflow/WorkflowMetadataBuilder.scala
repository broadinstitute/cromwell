package cromwell.engine.workflow

import akka.actor.ActorSystem
import cromwell.core.WorkflowId
import cromwell.database.SqlConverters._
import cromwell.database.obj.WorkflowExecution
import cromwell.engine
import cromwell.engine._
import cromwell.engine.backend.{OldStyleCallMetadata, OldStyleWorkflowDescriptor}
import cromwell.engine.db._
import cromwell.engine.finalcall.OldStyleFinalCall._
import cromwell.webservice._
import spray.json._
import wdl4s._

import scala.concurrent.Future

@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
object WorkflowMetadataBuilder {

  private type DBMap[+V] = Map[ExecutionDatabaseKey, V]

  private def build(workflowDescriptor: OldStyleWorkflowDescriptor,
                    execution: WorkflowExecution,
                    workflowOutputs: engine.WorkflowOutputs,
                    callMetadata: Map[FullyQualifiedName, Seq[OldStyleCallMetadata]]):
  WorkflowMetadataResponse = {

    val startDate = execution.startDt.toSystemOffsetDateTime
    val endDate = execution.endDt map { _.toSystemOffsetDateTime }
    val workflowInputs = workflowDescriptor.sourceFiles.inputsJson.parseJson.asInstanceOf[JsObject]

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
      calls = callMetadata)
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

@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
class WorkflowMetadataBuilder(id: WorkflowId, parameters: WorkflowMetadataQueryParameters)(implicit actorSystem: ActorSystem) {

  private[this] implicit val ec = actorSystem.dispatcher

  import WorkflowMetadataBuilder._
  import cromwell.engine.db.DataAccess.globalDataAccess

  private def futureAssertWorkflowExistsByState = globalDataAccess.assertWorkflowExistsByState(id)

  private def futureWorkflowExecutionAndAux = globalDataAccess.getWorkflowExecutionAndAux(id)

  // If outputs requested, don't retrieve the execution infos, only executions
  private def futureInfosByExecution = retrieveTraversable(parameters.outputs,
    globalDataAccess.infosByExecution(id),
    globalDataAccess.getExecutionsAsExecutionInfos(id))

  // If outputs requested, get the workflow outputs
  private def futureWorkflowOutputs = retrieveTraversable(parameters.outputs, globalDataAccess.getWorkflowOutputs(id))

  // If timings requested, we need the call statuses and execution events
  private def futureCallToStatusMap = retrieveMap(parameters.timings, globalDataAccess.getExecutionStatuses(id))

  // Drop the below if timings _or_ outputs requested
  private def futureCallInputs = retrieveTraversable(parameters.timings && parameters.outputs,
    globalDataAccess.getAllInputs(id))

  private def futureCallOutputs = retrieveTraversable(parameters.timings && parameters.outputs,
    globalDataAccess.getAllOutputs(id))

  private def futureCallCacheData = retrieveTraversable(parameters.timings && parameters.outputs,
    globalDataAccess.callCacheDataByExecution(id))

  private def futureRuntimeAttributes = retrieveMap(parameters.timings && parameters.outputs,
    globalDataAccess.getAllRuntimeAttributes(id))

  def build(): Future[WorkflowMetadataResponse] = {

    for {
      assertWorkflowExistsByState <- futureAssertWorkflowExistsByState
      workflowExecutionAndAux <- futureWorkflowExecutionAndAux
      infosByExecution <- futureInfosByExecution
      workflowOutputs <- futureWorkflowOutputs
      callToStatusMap <- futureCallToStatusMap
      callInputs <- futureCallInputs
      callOutputs <- futureCallOutputs
      callCacheData <- futureCallCacheData
      runtimeAttributes <- futureRuntimeAttributes

      // Database work complete, but we do need one more future to get the workflow descriptor
      workflowDescriptor <- workflowDescriptorFromExecutionAndAux(workflowExecutionAndAux)
    } yield {

      val execution = workflowExecutionAndAux.execution
      val nonFinalInfosByExecution = infosByExecution.filterNot(_.execution.callFqn.isFinalCall)

      val engineWorkflowOutputs = SymbolStoreEntry.toWorkflowOutputs(workflowOutputs)
      val callMetadata = CallMetadataBuilder.build(nonFinalInfosByExecution, callInputs, callOutputs,
        runtimeAttributes, callCacheData)

      WorkflowMetadataBuilder.build(workflowDescriptor, execution, engineWorkflowOutputs, callMetadata)
    }

  }
}
