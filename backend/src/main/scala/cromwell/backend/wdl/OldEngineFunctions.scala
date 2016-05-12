package cromwell.backend.wdl

import java.nio.file.Path

import cromwell.core.{OldCallContext, OldWorkflowContext}
import wdl4s.values._

import scala.language.postfixOps
import scala.util.{Success, Try}

/**
 * Default implementation of Wdl Standard Library functions executed at the workflow level.
 * This class is abstract to enforce that Backends extend it and customize the behavior if needed.
 */
@deprecated("Engine functions have now changed name and use Path contexts instead of String")
abstract class OldWorkflowEngineFunctions(val context: OldWorkflowContext) extends WdlStandardLibraryImpl {
  override def writeDirectory: Path = toPath(context.root)
}

/**
  * Provides a default implementation for Call specific engine functions that can be mixed in when implementing call engine function for a Backend.
  */
@deprecated("Engine functions have now changed name and use Path contexts instead of String")
trait OldCallEngineFunctions extends WdlStandardLibraryImpl {
  protected def stdout(context: OldCallContext): Try[WdlFile] = Success(WdlFile(context.stdout))
  protected def stderr(context: OldCallContext): Try[WdlFile] = Success(WdlFile(context.stderr))
}

