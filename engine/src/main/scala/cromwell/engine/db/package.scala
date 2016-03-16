package cromwell.engine

import cromwell.engine.ExecutionIndex.ExecutionIndex
import cromwell.engine.ExecutionStatus.ExecutionStatus
import cromwell.engine.backend.{BackendCallJobDescriptor, CallLogs}
import wdl4s.types.WdlFileType
import wdl4s.values.WdlFile

package object db {
  case class CallStatus(executionStatus: ExecutionStatus, returnCode: Option[Int], hash: Option[ExecutionHash], resultsClonedFrom: Option[BackendCallJobDescriptor]) {
    def isTerminal: Boolean = executionStatus.isTerminal
    def isStarting: Boolean = executionStatus == ExecutionStatus.Starting
  }

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

      var maybeStdout: Option[WdlFile] = None
      var maybeStderr: Option[WdlFile] = None
      var logs = Map.empty[String, WdlFile]

      entries foreach { entry =>
        val (keyFromSymbol, callLogName) = toCallLogKeyName(entry.key)
        // Skip all entries that don't belong to this key, have no .wdlValue, or are somehow the wrong type.
        if (keyFromSymbol == this && entry.wdlValue.isDefined && entry.wdlType == WdlFileType) {
          val wdlFile = entry.wdlValue.get.asInstanceOf[WdlFile]
          callLogName match {
            case StdoutSuffix => maybeStdout = Option(wdlFile)
            case StderrSuffix => maybeStderr = Option(wdlFile)
            case _ => logs += callLogName -> wdlFile
          }
        }
      }
      for {
        stdout <- maybeStdout
        stderr <- maybeStderr
      } yield CallLogs(stdout, stderr, Option(logs).filter(_.isEmpty))
    }
  }

  case class JesId(id: String)
  case class JesStatus(status: String)
}
