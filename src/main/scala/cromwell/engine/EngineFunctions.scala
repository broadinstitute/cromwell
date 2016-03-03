package cromwell.engine

import cromwell.engine.io.IoInterface
import wdl4s.values._

import scala.language.postfixOps
import scala.util.{Success, Try}

/**
 * Default implementation of Wdl Standard Library functions executed at the workflow level.
 * This class is abstract to enforce that Backends extend it and customize the behavior if needed.
 */
abstract class WorkflowEngineFunctions(override val interface: IoInterface, val context: WorkflowContext) extends WdlStandardLibraryImpl {
  override def tempFilePath: String = context.root
  override def globPath(glob: String): String = context.root
}

/**
  * Provides a default implementation for Call specific engine functions that can be mixed in when implementing call engine function for a Backend.
  */
trait CallEngineFunctions extends WdlStandardLibraryImpl {
  protected def stdout(context: CallContext): Try[WdlFile] = Success(WdlFile(context.stdout))
  protected def stderr(context: CallContext): Try[WdlFile] = Success(WdlFile(context.stderr))
}
