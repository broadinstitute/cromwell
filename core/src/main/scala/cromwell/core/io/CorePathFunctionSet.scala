package cromwell.core.io

import cromwell.core.CallContext
import cromwell.core.path.PathFactory
import cromwell.core.path.PathFactory.PathBuilders
import wom.expression.{IoFunctionSet, PathFunctionSet}

import scala.util.Try

class WorkflowCorePathFunctionSet(override val pathBuilders: PathBuilders) extends PathFunctionSet with PathFactory {
  private def fail(name: String) = throw new UnsupportedOperationException(s"$name is not implemented at the workflow level")
  override def sibling(of: String, path: String): String = buildPath(of).sibling(path).pathAsString
  override def isAbsolute(path: String): Boolean = Try(buildPath(path)).map(_.isAbsolute).toOption.contains(true)
  override def name(path: String) = buildPath(path).name

  // Call level functions
  override def relativeToHostCallRoot(path: String): String = fail("relativeToHostCallRoot")
  override def stdout: String = fail("stdout")
  override def stderr: String = fail("stderr")
}

class CallCorePathFunctionSet(pathBuilders: PathBuilders, callContext: CallContext) extends WorkflowCorePathFunctionSet(pathBuilders) {
  override def relativeToHostCallRoot(path: String) = if (isAbsolute(path)) path else callContext.root.resolve(path).pathAsString
  override def stdout = callContext.standardPaths.output.pathAsString
  override def stderr = callContext.standardPaths.error.pathAsString
}

trait WorkflowCorePathFunctions extends { this: IoFunctionSet with PathFactory =>
  override lazy val pathFunctions = new WorkflowCorePathFunctionSet(pathBuilders)
}


trait CallCorePathFunctions extends { this: IoFunctionSet with PathFactory =>
  def callContext: CallContext
  override lazy val pathFunctions = new CallCorePathFunctionSet(pathBuilders, callContext)
}
