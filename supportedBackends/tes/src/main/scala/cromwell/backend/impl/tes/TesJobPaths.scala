package cromwell.backend.impl.tes

import java.nio.file.{Path, Paths}

import com.typesafe.config.Config
import cromwell.backend.io.{JobPaths, WorkflowPaths}
import cromwell.backend.{BackendJobDescriptorKey, BackendWorkflowDescriptor}
import cromwell.core.path.PathBuilder
import wdl4s.values.{WdlArray, WdlFile, WdlMap, WdlValue}

class TesJobPaths(val jobKey: BackendJobDescriptorKey,
                  workflowDescriptor: BackendWorkflowDescriptor,
                  config: Config,
                  pathBuilders: List[PathBuilder] = WorkflowPaths.DefaultPathBuilders) extends TesWorkflowPaths(
  workflowDescriptor, config, pathBuilders) with JobPaths {
  import JobPaths._

  override lazy val callExecutionRoot = { callRoot.resolve("execution") }
  val callDockerRoot = callPathBuilder(dockerWorkflowRoot, jobKey)
  val callExecutionDockerRoot = callDockerRoot.resolve("execution")
  val callInputsRoot = callRoot.resolve("inputs")
  var containerWorkingDir = callExecutionDockerRoot

  // Utility for converting a WdlValue so that the path is localized to the
  // container's filesystem.
  def toContainerPath(path: WdlValue): WdlValue = {
    path match {
      case file: WdlFile => {
        val localPath = Paths.get(file.valueString).toAbsolutePath
        val containerPath = containerInput(localPath.toString)
        WdlFile(containerPath)
      }
      case array: WdlArray => WdlArray(array.wdlType, array.value map toContainerPath)
      case map: WdlMap => WdlMap(map.wdlType, map.value mapValues toContainerPath)
      case wdlValue => wdlValue
    }
  }

  private def prefixScheme(path: String): String = "file://" + path

  def storageInput(path: String): String = prefixScheme(path)

  // Given an output path, return a path localized to the storage file system
  def storageOutput(path: String): String = {
    prefixScheme(callExecutionRoot.resolve(path).toString)
  }

  def containerInput(path: String): String = {
    callDockerRoot.resolve("inputs").toString + cleanPathForContainer(Paths.get(path))
  }

  // Given an output path, return a path localized to the container file system
  def containerOutput(path: String): String = containerExec(path)

  // TODO this could be used to create a separate directory for outputs e.g.
  // callDockerRoot.resolve("outputs").resolve(name).toString

  // Given an file name, return a path localized to the container's execution directory
  def containerExec(name: String): String = {
    containerWorkingDir.resolve(cleanPathForContainer(Paths.get(name))).toString
  }

  def cleanPathForContainer(path: Path): String = {
    path.toAbsolutePath match {
      case p if p.startsWith(executionRoot) => {
        /* For example:
          *
          * p = /abs/path/to/cromwell-executions/three-step/f00ba4/call-ps/stdout.txt
          * localExecutionRoot = /abs/path/to/cromwell-executions
          * subpath = three-step/f00ba4/call-ps/stdout.txt
          *
          * return value = /root/three-step/f00ba4/call-ps/stdout.txt
          *
          */
        val subpath = p.subpath(executionRoot.getNameCount, p.getNameCount)
        subpath.toString
      }
      case _ => path.toString
    }
  }
}
