package cromwell.engine.db

import cromwell.database.SqlConverters._
import cromwell.database.obj.{Execution, Symbol}
import cromwell.engine.ExecutionIndex._
import cromwell.engine.workflow.BackendCallKey
import cromwell.engine.{ExecutionStatus, SymbolStoreEntry, SymbolStoreKey}
import wdl4s.types.{WdlPrimitiveType, WdlType}
import wdl4s.values.{SymbolHash, WdlValue}
import wdl4s.{NamespaceWithWorkflow, ReportableSymbol, Scatter}

import scala.util.{Failure, Success, Try}

object EngineConverters {

  private def wdlValueToDbValue(v: WdlValue): String = v.wdlType match {
    case p: WdlPrimitiveType => v.valueString
    case o => v.toWdlString
  }

  private def dbEntryToWdlValue(dbValue: String, wdlType: WdlType): WdlValue = wdlType match {
    // .get here is because we trust the value in the database is coercible to the given type
    case p: WdlPrimitiveType => p.coerceRawValue(dbValue).get
    case o => wdlType.fromWdlString(dbValue)
  }

  val IoInput = "INPUT"
  val IoOutput = "OUTPUT"

  implicit class SymbolToSymbolStoreEntry(val symbol: Symbol) extends AnyVal {
    def toSymbolStoreEntry: SymbolStoreEntry = {
      val wdlType = WdlType.fromWdlString(symbol.wdlType)
      val value = symbol.wdlValue map { v =>
        dbEntryToWdlValue(v.toRawString, wdlType)
      } orElse Option(dbEntryToWdlValue("", wdlType))

      new SymbolStoreEntry(
        new SymbolStoreKey(
          symbol.scope,
          symbol.name,
          symbol.index.toIndex,
          input = symbol.io == IoInput // input = true, if db contains "INPUT"
        ),
        wdlType,
        value,
        symbol.symbolHash map SymbolHash
      )
    }
  }

  implicit class SymbolStoreEntryToSymbol(val symbolStoreEntry: SymbolStoreEntry) extends AnyVal {
    def toSymbol(workflowOutputs: Seq[ReportableSymbol])(workflowExecutionId: Int) = {
      val reportableResult = workflowOutputs exists {
        _.fullyQualifiedName == symbolStoreEntry.key.fqn
      }
      val value = symbolStoreEntry.wdlValue map wdlValueToDbValue flatMap {
        _.toNonEmptyClob
      }
      new Symbol(
        workflowExecutionId,
        symbolStoreEntry.key.scope,
        symbolStoreEntry.key.name,
        symbolStoreEntry.key.index.fromIndex,
        if (symbolStoreEntry.key.input) IoInput else IoOutput,
        reportableResult,
        symbolStoreEntry.wdlType.toWdlString,
        value,
        symbolStoreEntry.symbolHash map {
          _.value
        }
      )
    }
  }

  implicit class EnhancedExecution(val execution: Execution) extends AnyVal {
    def isShard = execution.index.toIndex.isShard

    def isScatter = execution.callFqn.contains(Scatter.FQNIdentifier)

    def isCollector(keys: Traversable[Execution]): Boolean = {
      !isShard && (keys exists { e => (e.execution.callFqn == execution.callFqn) && e.isShard })
    }

    def toKey = ExecutionDatabaseKey(execution.callFqn, execution.index.toIndex, execution.attempt)

    def executionStatus = ExecutionStatus.withName(execution.status)

    def toBackendCallKey(ns: NamespaceWithWorkflow): Try[BackendCallKey] = {
      val call = ns.workflow.calls.find(_.fullyQualifiedName == execution.callFqn)
      call.map(BackendCallKey(_, execution.index.toIndex, execution.attempt)) match {
        case Some(key) => Success(key)
        case None => Failure(new IllegalStateException(s"Could not find call with FQN '$execution.callFqn'"))
      }
    }
  }

}
