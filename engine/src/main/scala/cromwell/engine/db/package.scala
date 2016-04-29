package cromwell.engine

import cromwell.engine.ExecutionIndex.ExecutionIndex
import cromwell.engine.ExecutionStatus.ExecutionStatus
import cromwell.engine.backend.{OldStyleBackendCallJobDescriptor, ExecutionHash}
import wdl4s.Scatter

package object db {
  @deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
  case class CallStatus(executionStatus: ExecutionStatus, returnCode: Option[Int], hash: Option[ExecutionHash], resultsClonedFrom: Option[OldStyleBackendCallJobDescriptor]) {
    def isTerminal: Boolean = executionStatus.isTerminal
    def isStarting: Boolean = executionStatus == ExecutionStatus.Starting
  }

  // Uniquely identify an entry in the execution table
  @deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
  case class ExecutionDatabaseKey(fqn: FullyQualifiedName, index: ExecutionIndex, attempt: Int) {
    def isCollector(keys: Traversable[ExecutionDatabaseKey]): Boolean = {
      index.isEmpty &&
        (keys exists { e =>
          (e.fqn == fqn) && e.index.isDefined
        })
    }

    def isScatter = fqn.contains(Scatter.FQNIdentifier)
  }

  @deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
  case class JesId(id: String)
  @deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
  case class JesStatus(status: String)
}
