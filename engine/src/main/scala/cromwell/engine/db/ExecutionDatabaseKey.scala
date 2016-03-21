package cromwell.engine.db

import cromwell.engine.ExecutionIndex._
import cromwell.engine._
import cromwell.engine.backend._
import wdl4s.types.WdlFileType
import wdl4s.values.WdlFile

import scala.collection.immutable.Seq

object ExecutionDatabaseKey {
  val CallLogPrefix = "$log"
  val StdoutSuffix = "stdout"
  val StderrSuffix = "stderr"

  def isCallLog(entry: SymbolStoreEntry): Boolean = isCallLog(entry.key)

  def isCallLog(key: SymbolStoreKey): Boolean = key.name.startsWith(CallLogPrefix + "_")

  def toCallLogKey(entry: SymbolStoreEntry): ExecutionDatabaseKey = toCallLogKey(entry.key)

  def toCallLogKey(key: SymbolStoreKey): ExecutionDatabaseKey = toCallLogKeyName(key)._1

  def toCallLogKeyName(key: SymbolStoreKey): (ExecutionDatabaseKey, String) = {
    key.name.split("_", 3) match {
      case Array(CallLogPrefix, attemptString, callLogName) =>
        (ExecutionDatabaseKey(key.scope, key.index, attemptString.toInt), callLogName)
      case unexpected =>
        throw new Exception(s"Unexpected call log key ('${unexpected.mkString("','")}')")
    }
  }

  def toAttemptedCallLogs(maybeLogs: Seq[(ExecutionDatabaseKey, Option[CallLogs])]): AttemptedCallLogs = {
    val logs = maybeLogs collect { case (key, Some(callLog)) => key -> callLog }
    val groupByIndex = logs groupBy { case (key, _) => key.index }
    val sortedByIndex = groupByIndex.toIndexedSeq sortBy { case (index, _) => index }
    val attemptedCallLogs = sortedByIndex map { case (_, callLogs) => mapSortByAttempt(callLogs) }
    attemptedCallLogs
  }

  private def mapSortByAttempt(logs: Seq[(ExecutionDatabaseKey, CallLogs)]): Seq[CallLogs] = {
    val sortedByAttempt = logs sortBy { case (key, callLog) => key.attempt }
    val mappedToLogs = sortedByAttempt map { case (_, callLog) => callLog }
    mappedToLogs
  }
}

// Uniquely identify an entry in the execution table
case class ExecutionDatabaseKey(fqn: FullyQualifiedName, index: ExecutionIndex, attempt: Int) {
  def isCollector(keys: Traversable[ExecutionDatabaseKey]): Boolean = {
    index.isEmpty &&
      (keys exists { e =>
        (e.fqn == fqn) && e.index.isDefined
      })
  }

  def toCallOutputs(callLogs: CallLogs): CallOutputs = {
    import ExecutionDatabaseKey._
    val callLogMap = Map(StdoutSuffix -> callLogs.stdout, StderrSuffix -> callLogs.stderr) ++
      callLogs.backendLogs.getOrElse(Map.empty)

    callLogMap map {
      case (suffix, value) => s"${CallLogPrefix}_${attempt}_$suffix" -> CallOutput(value, None)
    }
  }

  def toCallLogs(entries: Traversable[SymbolStoreEntry]): Option[CallLogs] = {
    import ExecutionDatabaseKey._

    case class Accumulator(stdout: Option[WdlFile], stderr: Option[WdlFile], logs: Map[String, WdlFile])

    val acc = entries.foldLeft(Accumulator(None, None, Map.empty)) { (acc, entry) =>
      val (keyFromSymbol, callLogName) = toCallLogKeyName(entry.key)
      // Skip all entries that don't belong to this key, have no .wdlValue, or are somehow the wrong type.
      if (keyFromSymbol == this && entry.wdlValue.isDefined && entry.wdlType == WdlFileType) {
        val wdlFile = entry.wdlValue.get.asInstanceOf[WdlFile]
        callLogName match {
          case StdoutSuffix => acc.copy(stdout = Option(wdlFile))
          case StderrSuffix => acc.copy(stderr = Option(wdlFile))
          case _ => acc.copy(logs = acc.logs + (callLogName -> wdlFile))
        }
      } else acc
    }

    for {
      stdout <- acc.stdout
      stderr <- acc.stderr
    } yield CallLogs(stdout, stderr, Option(acc.logs).filterNot(_.isEmpty))
  }
}
