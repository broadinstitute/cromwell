package cromwell.backend.impl.tes

import java.nio.file.{FileSystems, Paths}

import cromwell.backend.io.JobPaths
import cromwell.backend.sfs.SharedFileSystemExpressionFunctions
import cromwell.backend.{BackendConfigurationDescriptor, BackendJobDescriptor}
import wdl4s.parser.MemoryUnit
import wdl4s.util.TryUtil
import wdl4s.values.{WdlArray, WdlFile, WdlMap, WdlSingleFile, WdlValue}

import scala.util.Try

final case class TesTaskMessage(
                          name: String,
                          projectId: String,
                          description: Option[String],
                          inputs: Option[Seq[TaskParameter]],
                          outputs: Option[Seq[TaskParameter]],
                          resources: Resources,
                          taskId: String,
                          docker: Seq[DockerExecutor]
                        )

final case class DockerExecutor(
                                 imageName: String,
                                 cmd: Seq[String],
                                 workDir: Option[String],
                                 stdout: String,
                                 stderr: String,
                                 stdin: Option[String]
                               )

final case class TaskParameter(
                                name: String,
                                description: Option[String],
                                location: String,
                                path: String,
                                `class`: String,
                                create: Boolean
                              )

final case class Resources(
                            minimumCpuCores: Option[Int],
                            preemtible: Option[Boolean],
                            minimumRamGb: Option[Int],
                            volumes: Seq[Volume],
                            zones: Option[Seq[String]]
                          )

final case class Volume(
                         name: Option[String],
                         sizeGb: Option[Int],
                         source: Option[String],
                         mountPoint: Option[String]
                       )

final case class TesTask(jobDescriptor: BackendJobDescriptor,
                         configurationDescriptor: BackendConfigurationDescriptor) {

  private val workflowDescriptor = jobDescriptor.workflowDescriptor
  private val jobPaths = new JobPaths(workflowDescriptor, configurationDescriptor.backendConfig, jobDescriptor.key)
  private val callEngineFunction = SharedFileSystemExpressionFunctions(jobPaths, List(FileSystems.getDefault))

  val name = jobDescriptor.call.task.name
  val taskId = jobDescriptor.toString
  val project = jobDescriptor.workflowDescriptor.workflowOptions.getOrElse(
    "project",
    jobDescriptor.call.rootWorkflow.unqualifiedName
  )
  val desc = None

  private val runtimeAttributes = {
    val lookup = jobDescriptor.inputs.apply _
    val evaluateAttrs = jobDescriptor.call.task.runtimeAttributes.attrs mapValues (_.evaluate(lookup, callEngineFunction))
    // Fail the call if runtime attributes can't be evaluated
    val runtimeMap = TryUtil.sequenceMap(evaluateAttrs, "Runtime attributes evaluation").get
    TesRuntimeAttributes(runtimeMap, jobDescriptor.workflowDescriptor.workflowOptions)
  }

  private def toDockerPath(path: WdlValue): WdlValue = {
    path match {
      case file: WdlFile => {
        val localPath = Paths.get(file.valueString)

        localPath.toAbsolutePath match {
          case p if p.startsWith(jobPaths.DockerRoot) => WdlFile(p.toString)
          case p =>
            val fileName = p.getFileName
            val callPath = jobPaths.callRoot.resolve(fileName)
            val subPath = callPath.subpath(jobPaths.executionRoot.getNameCount, callPath.getNameCount)
            WdlFile(jobPaths.DockerRoot.resolve(subPath).toString)
        }
      }
      case array: WdlArray => WdlArray(array.wdlType, array.value map toDockerPath)
      case map: WdlMap => WdlMap(map.wdlType, map.value mapValues toDockerPath)
      case wdlValue => wdlValue
    }
  }

  val inputs: Seq[TaskParameter] = {
    jobDescriptor.inputs.toSeq.map {
      case (varName, f: WdlSingleFile) => TaskParameter(
        varName,
        Option("description"),
        f.value,
        toDockerPath(f).toString,
        "file",
        false
      )
    }
  }

  val outputs = Seq()

  //  private def getCommonVolumes(inputs: Seq[Map[String,String]],
  //                               outputs: Seq[Map[String,String]]): Seq[Map[String, String]] = {
  //
  //  }

  val command: Try[String] = {
    jobDescriptor
      .call
      .instantiateCommandLine(
        jobDescriptor.inputs,
        callEngineFunction,
        toDockerPath
      )
  }

  // TODO - command shouldn't be wrapped in a sub-shell like this
  val dockerExecutor = for {
    cmd <- command
  } yield DockerExecutor(
    runtimeAttributes.dockerImage.get,
    Seq("/bin/bash", "-c", cmd),
    runtimeAttributes.dockerWorkingDir,
    "/tmp/stdout",
    "/tmp/stderr",
    None
  )

  // TODO - resolve TES schema around memory format Int -> Double
  val resources = Resources(
    Some(runtimeAttributes.cpu),
    None,
    Some(runtimeAttributes.memory.to(MemoryUnit.GB).amount.toInt),
    Seq(
      Volume(
        Some("cromwell_inputs"),
        Some(runtimeAttributes.disk.to(MemoryUnit.GB).amount.toInt),
        None,
        Some("/tmp")
      )
    ),
    None
  )
}
