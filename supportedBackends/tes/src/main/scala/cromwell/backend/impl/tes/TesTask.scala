package cromwell.backend.impl.tes

import java.nio.file.{FileSystems, Paths}

import cromwell.backend.io.JobPaths
import cromwell.backend.sfs.SharedFileSystemExpressionFunctions
import cromwell.backend.{BackendConfigurationDescriptor, BackendJobDescriptor}
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
                            minimumRamGb: Option[Double],
                            volumes: Seq[Volume],
                            zones: Option[Seq[String]]
                          )

final case class Volume(
                         name: Option[String],
                         sizeGb: Option[Double],
                         source: Option[String],
                         mountPoint: Option[String]
                       )

final case class TesTask(jobDescriptor: BackendJobDescriptor,
                         configurationDescriptor: BackendConfigurationDescriptor) {

  val workflowDescriptor = jobDescriptor.workflowDescriptor
  val jobPaths = new JobPaths(workflowDescriptor, configurationDescriptor.backendConfig, jobDescriptor.key)
  val callEngineFunction = SharedFileSystemExpressionFunctions(jobPaths, List(FileSystems.getDefault))

  val name = jobDescriptor.call.task.name
  val taskId = jobDescriptor.toString
  val project = jobDescriptor.workflowDescriptor.workflowOptions.getOrElse(
    "project",
    jobDescriptor.call.rootWorkflow.unqualifiedName
  )
  val desc = None

  val runtimeAttributes = {
    val lookup = jobDescriptor.inputs.apply _
    val evaluateAttrs = jobDescriptor.call.task.runtimeAttributes.attrs mapValues (_.evaluate(lookup, callEngineFunction))
    // Fail the call if runtime attributes can't be evaluated
    val runtimeMap = TryUtil.sequenceMap(evaluateAttrs, "Runtime attributes evaluation").get
    TesRuntimeAttributes(runtimeMap, jobDescriptor.workflowDescriptor.workflowOptions)
  }

  private def toDockerPath(path: WdlValue): WdlValue = {
    path match {
      case file: WdlFile => WdlFile(jobPaths.toDockerPath(Paths.get(path.valueString)).toString)
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
  //  def prepareOutputs(jobDescriptor: BackendJobDescriptor): Seq[Map[String, String]] = {
  //    val wdlFileOutputs = jobDescriptor.key.scope.task.outputs flatMap { taskOutput =>
  //      taskOutput.requiredExpression.evaluateFiles(lookup, NoFunctions, taskOutput.wdlType) match {
  //        case Success(wdlFiles) => wdlFiles
  //        case Failure(ex) =>
  //          jobLogger.warn(s"Could not evaluate $taskOutput: ${ex.getMessage}", ex)
  //          Seq.empty[WdlFile]
  //      }
  //    }
  //
  //    wdlFileOutputs.distinct map { wdlFile =>
  //      val lPath = wdlFile match {
  //        case WdlSingleFile(filePath) => jobPaths.callExecutionRoot.resolve(filePath).toString
  //        case WdlGlobFile(filePath) => jobPaths.callExecutionRoot.resolve(filePath).toString
  //      }
  //      val mPath = toDockerPath(wdlFile).toString
  //      Map("name" -> "",
  //          "localPath" -> lPath,
  //          "mountPath" -> mPath)
  //    }
  //  }

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
    Option(runtimeAttributes.cpu),
    None,
    Option(runtimeAttributes.memory.amount.toInt),
    Seq(
      Volume(
        Some("cromwell_inputs"),
        Some(runtimeAttributes.disk.amount.toInt),
        None,
        Some("/tmp")
      )
    ),
    None
  )
}
