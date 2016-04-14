package cromwell.engine.db

import cromwell.engine._
import cromwell.engine.backend._
import cromwell.engine.db.slick.{Execution, ExecutionInfo}
import cromwell.engine.finalcall.FinalCall
import wdl4s.values.WdlFile

case class ExecutionInfosByExecution(execution: Execution, executionInfos: Seq[ExecutionInfo]) {
  import ExecutionInfosByExecution._
  lazy val callLogs: Option[CallLogs] = {
    executionInfos.foldLeft(Accumulator(None, None, Map.empty))(accumulateExecutionInfo) match {
      case Accumulator(Some(stdout), None, _) =>
        throw new IllegalArgumentException(s"stderr was not found for stdout $stdout")
      case Accumulator(None, Some(stderr), _) =>
        throw new IllegalArgumentException(s"stdout was not found for stderr $stderr")
      case Accumulator(None, None, logs) if logs.nonEmpty =>
        throw new IllegalArgumentException(s"stdout and stderr were empty logs $logs")
      case acc =>
        for {
          stdout <- acc.stdout
          stderr <- acc.stderr
        } yield CallLogs(stdout, stderr, Option(acc.logs).filterNot(_.isEmpty))
    }
  }
}

object ExecutionInfosByExecution {
  private val CallLogPrefix = "$log"
  private val StdoutSuffix = "stdout"
  private val StderrSuffix = "stderr"

  private case class Accumulator(stdout: Option[WdlFile], stderr: Option[WdlFile], logs: Map[String, WdlFile])

  private def accumulateExecutionInfo(acc: Accumulator, executionInfo: ExecutionInfo): Accumulator = {
    executionInfo.key.split("_", 2) match {
      case Array(CallLogPrefix, StdoutSuffix) => acc.copy(stdout = executionInfo.value map { WdlFile(_) })
      case Array(CallLogPrefix, StderrSuffix) => acc.copy(stderr = executionInfo.value map { WdlFile(_) })
      case Array(CallLogPrefix, callLogName) if executionInfo.value.isDefined =>
        acc.copy(logs = acc.logs + (callLogName -> WdlFile(executionInfo.value.get)))
      case _ => acc
    }
  }

  def toCallLogMap(callLogs: CallLogs): Map[String, Option[String]] = {
    val callLogMap = Map(StdoutSuffix -> callLogs.stdout, StderrSuffix -> callLogs.stderr) ++
      callLogs.backendLogs.getOrElse(Map.empty)

    callLogMap map {
      case (suffix, value) => s"${CallLogPrefix}_$suffix" -> Option(value.valueString)
    }
  }

  def toWorkflowLogs(executionInfosByExecutions: Traversable[ExecutionInfosByExecution]): WorkflowLogs = {
    import FinalCall._
    executionInfosByExecutions filter {
      !_.execution.toKey.isFinalCall
    } groupBy {
      _.execution.callFqn
    } mapValues {
      toAttemptedCallLogs
    } filterNot {
      case (_, attemptedCallLogs) => attemptedCallLogs.isEmpty
    }
  }

  private def toAttemptedCallLogs(executionInfosByExecutions: Traversable[ExecutionInfosByExecution]):
  AttemptedCallLogs = {
    toAttemptedCallLogs(executionInfosByExecutions.toIndexedSeq map toKeyCallLogs)
  }

  private def toKeyCallLogs(executionInfosByExecution: ExecutionInfosByExecution):
  (ExecutionDatabaseKey, Option[CallLogs]) = {
    import ExecutionIndex._
    val execution = executionInfosByExecution.execution
    val key = ExecutionDatabaseKey(execution.callFqn, execution.index.toIndex, execution.attempt)
    val callLogs = executionInfosByExecution.callLogs
    key -> callLogs
  }

  private def toAttemptedCallLogs(maybeLogs: Seq[(ExecutionDatabaseKey, Option[CallLogs])]): AttemptedCallLogs = {
    val logs = maybeLogs collect { case (key, Some(callLog)) => key -> callLog }
    val groupByIndex = logs groupBy { case (key, _) => key.index }
    val sortedByIndex = groupByIndex.toIndexedSeq sortBy { case (index, _) => index }
    val attemptedCallLogs = sortedByIndex map { case (_, callLogs) => mapSortByAttempt(callLogs) }
    attemptedCallLogs map { _.toIndexedSeq }
  }

  private def mapSortByAttempt(logs: Seq[(ExecutionDatabaseKey, CallLogs)]): Seq[CallLogs] = {
    val sortedByAttempt = logs sortBy { case (key, callLog) => key.attempt }
    val mappedToLogs = sortedByAttempt map { case (_, callLog) => callLog }
    mappedToLogs
  }

  /**
    * Group by execution, then remove the execution keys of the execution -> execution info tuple.
    * The net result is Execution -> Seq[ExecutionInfo].
    */
  def fromRawTuples(rawTuples: Seq[(Execution, ExecutionInfo)]): Seq[ExecutionInfosByExecution] = {
    val groupedTuples = rawTuples groupBy {
      case (execution, _) => execution
    } mapValues {
      _ map { case (_, executionInfo) => executionInfo }
    }
    val infosByExecution = groupedTuples map (ExecutionInfosByExecution.apply _).tupled
    infosByExecution.toSeq
  }
}
