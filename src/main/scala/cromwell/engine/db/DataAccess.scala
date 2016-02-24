package cromwell.engine.db

import cromwell.engine.ExecutionStatus.ExecutionStatus
import cromwell.engine.backend.CallLogs
import cromwell.engine.db.slick._
import cromwell.engine.workflow.{BackendCallKey, ExecutionStoreKey}
import cromwell.engine.{WorkflowOutputs, _}
import cromwell.webservice.{CallCachingParameters, WorkflowQueryParameters, WorkflowQueryResponse}
import wdl4s.values.WdlFile
import wdl4s.{CallInputs, _}

import scala.concurrent.{ExecutionContext, Future}
import scala.language.postfixOps

object DataAccess {
  val globalDataAccess: DataAccess = new slick.SlickDataAccess()
  case class ExecutionKeyToJobKey(executionKey: ExecutionDatabaseKey, jobKey: JobKey)

  val IsStandardStream: SymbolStoreEntry => Boolean = s => s.key.fqn.endsWith(".$stdout") || s.key.fqn.endsWith(".$stderr")

  val StdoutSuffix = ".$stdout"
  val StderrSuffix = ".$stderr"
}

trait DataAccess extends AutoCloseable {
  import DataAccess._

  def createWorkflow(workflowDescriptor: WorkflowDescriptor,
                     workflowInputs: Traversable[SymbolStoreEntry],
                     calls: Traversable[Scope])(implicit ec: ExecutionContext): Future[Unit]

  def getWorkflowState(workflowId: WorkflowId)(implicit ec: ExecutionContext): Future[Option[WorkflowState]]

  def getWorkflow(workflowId: WorkflowId)(implicit ec: ExecutionContext): Future[WorkflowDescriptor]

  def getWorkflow(workflowExecutionId: Int)(implicit ec: ExecutionContext): Future[WorkflowDescriptor]

  def getWorkflowsByState(states: Traversable[WorkflowState])
                         (implicit ec: ExecutionContext): Future[Traversable[WorkflowDescriptor]]

  def getExecutionInfos(workflowId: WorkflowId, call: Call, attempt: Int)(implicit ec: ExecutionContext): Future[Traversable[ExecutionInfo]]

  def updateExecutionInfo(workflowId: WorkflowId, callKey: BackendCallKey, key: String, value: Option[String])
                         (implicit ec: ExecutionContext): Future[Unit]

  def updateWorkflowState(workflowId: WorkflowId, workflowState: WorkflowState)
                         (implicit ec: ExecutionContext): Future[Unit]

  def getAllSymbolStoreEntries(workflowId: WorkflowId)
                              (implicit ec: ExecutionContext): Future[Traversable[SymbolStoreEntry]]

  // TODO needed to support compatibility with current code, this seems like an inefficient way of getting
  // TODO workflow outputs.
  /** Returns all outputs for this workflowId */
  def getWorkflowOutputs(workflowId: WorkflowId)(implicit ec: ExecutionContext): Future[Traversable[SymbolStoreEntry]]

  protected def getAllOutputs(workflowId: WorkflowId)(implicit ec: ExecutionContext): Future[Traversable[SymbolStoreEntry]]

  def getOutputs(workflowId: WorkflowId)(implicit ec: ExecutionContext): Future[Traversable[SymbolStoreEntry]] = {
    getAllOutputs(workflowId)(ec) map { _ filterNot IsStandardStream }
  }

  def getAllInputs(workflowId: WorkflowId)(implicit ec: ExecutionContext): Future[Traversable[SymbolStoreEntry]]

  /** Get all output symbols for the scope of this call including standard stream symbols. */
  protected def getAllOutputs(workflowId: WorkflowId, key: ExecutionDatabaseKey)
                             (implicit ec: ExecutionContext): Future[Traversable[SymbolStoreEntry]]

  /** Get all "real" outputs for the scope of this call (i.e. exclude standard stream symbols). */
  def getOutputs(workflowId: WorkflowId, key: ExecutionDatabaseKey)
                (implicit ec: ExecutionContext): Future[Traversable[SymbolStoreEntry]] = {
    getAllOutputs(workflowId, key)(ec) map { _ filterNot IsStandardStream }
  }

  /** Get the standard stream symbols for the scope of this call. */
  def getStandardStreamSymbols(workflowId: WorkflowId, key: ExecutionDatabaseKey)
                              (implicit ec: ExecutionContext): Future[CallLogs] = {

    // Folds stdout and stderr symbols into an accumulator CallLogs object.
    def foldLogs(acc: CallLogs, symbol: SymbolStoreEntry): CallLogs = {
      symbol match {
        case s if s.key.name.endsWith(StdoutSuffix) => acc.copy(stdout = s.wdlValue.get.asInstanceOf[WdlFile])
        case s if s.key.name.endsWith(StderrSuffix) => acc.copy(stderr = s.wdlValue.get.asInstanceOf[WdlFile])
        case s => acc
      }
    }

    // Initial "zero" value of a CallLogs with blank stdout and stderr.
    val zeroLogs = CallLogs(stdout = WdlFile(""), stderr = WdlFile(""))

    getAllOutputs(workflowId, key)(ec) map { symbols =>
      symbols.foldLeft(zeroLogs) { case (acc, symbol) => foldLogs(acc, symbol) }
    }
  }

  /** Get all inputs for the scope of this call. */
  def getInputs(id: WorkflowId, call: Call)(implicit ec: ExecutionContext): Future[Traversable[SymbolStoreEntry]]

  /** Should fail if a value is already set.  The keys in the Map are locally qualified names. */
  def setOutputs(workflowId: WorkflowId, key: ExecutionDatabaseKey, callOutputs: WorkflowOutputs,
                 workflowOutputFqns: Seq[ReportableSymbol])(implicit ec: ExecutionContext): Future[Unit]

  /** Updates the existing input symbols to replace expressions with real values **/
  def updateCallInputs(workflowId: WorkflowId, key: BackendCallKey, callInputs: CallInputs)
                      (implicit ec: ExecutionContext): Future[Int]

  def setExecutionEvents(workflowId: WorkflowId, callFqn: String, shardIndex: Option[Int], attempt: Int,
                         events: Seq[ExecutionEventEntry])(implicit ec: ExecutionContext): Future[Unit]

  /** Gets a mapping from call FQN to an execution event entry list */
  def getAllExecutionEvents(workflowId: WorkflowId)
                           (implicit ec: ExecutionContext): Future[Map[ExecutionDatabaseKey, Seq[ExecutionEventEntry]]]

  def setStatus(workflowId: WorkflowId, keys: Traversable[ExecutionDatabaseKey], executionStatus: ExecutionStatus)
               (implicit ec: ExecutionContext): Future[Unit] = {
    setStatus(workflowId, keys, CallStatus(executionStatus, None, None, None))
  }

  def setStatus(workflowId: WorkflowId, keys: Traversable[ExecutionDatabaseKey], callStatus: CallStatus)
               (implicit ec: ExecutionContext): Future[Unit]

  def getExecutionStatuses(workflowId: WorkflowId)
                          (implicit ec: ExecutionContext): Future[Map[ExecutionDatabaseKey, CallStatus]]

  /** Return all execution entries for the FQN, including collector and shards if any */
  def getExecutionStatuses(workflowId: WorkflowId, fqn: FullyQualifiedName)
                          (implicit ec: ExecutionContext): Future[Map[ExecutionDatabaseKey, CallStatus]]

  def getExecutionStatus(workflowId: WorkflowId, key: ExecutionDatabaseKey)
                        (implicit ec: ExecutionContext): Future[Option[CallStatus]]

  def insertCalls(workflowId: WorkflowId, keys: Traversable[ExecutionStoreKey])
                 (implicit ec: ExecutionContext): Future[Unit]

  def getExecutions(id: WorkflowId)(implicit ec: ExecutionContext): Future[Traversable[Execution]]

  def getExecutionsForRestart(id: WorkflowId)(implicit ec: ExecutionContext): Future[Traversable[Execution]]

  def getExecutionsWithResuableResultsByHash(hash: String)
                                            (implicit ec: ExecutionContext): Future[Traversable[Execution]]

  /** Fetch the workflow having the specified `WorkflowId`. */
  def getWorkflowExecution(workflowId: WorkflowId)(implicit ec: ExecutionContext): Future[WorkflowExecution]

  def getWorkflowExecutionAux(id: WorkflowId)(implicit ec: ExecutionContext): Future[WorkflowExecutionAux]

  def updateWorkflowOptions(workflowId: WorkflowId, workflowOptionsJson: String)(implicit ec: ExecutionContext): Future[Unit]

  def resetNonResumableExecutions(workflowId: WorkflowId, isResumable: ExecutionAndExecutionInfo => Boolean)(implicit ec: ExecutionContext): Future[Unit]

  def findResumableExecutions(workflowId: WorkflowId,
                              isResumable: ExecutionAndExecutionInfo => Boolean,
                              jobKeyBuilder: ExecutionAndExecutionInfo => JobKey)
                             (implicit ec: ExecutionContext): Future[Traversable[ExecutionKeyToJobKey]]

  def queryWorkflows(queryParameters: WorkflowQueryParameters)
                    (implicit ec: ExecutionContext): Future[WorkflowQueryResponse]

  def updateCallCaching(cachingParameters: CallCachingParameters)(implicit ec: ExecutionContext): Future[Int]

  def infosByExecution(id: WorkflowId)(implicit ec: ExecutionContext): Future[Traversable[ExecutionInfosByExecution]]

  def writeStandardStreamSymbols(id: WorkflowId, key: ExecutionDatabaseKey, stdout: String, stderr: String)(implicit ec: ExecutionContext): Future[Unit] = {
    // Write the magical $stdout and $stderr symbols.
    val outputs = Map(StdoutSuffix -> stdout, StderrSuffix -> stderr) map {
      case (suffix, value) => key.fqn + suffix -> CallOutput(WdlFile(value), None) }
    setOutputs(id, key, outputs, Seq.empty)
  }
}
