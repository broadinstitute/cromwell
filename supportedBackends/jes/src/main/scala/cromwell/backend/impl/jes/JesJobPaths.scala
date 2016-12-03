package cromwell.backend.impl.jes

import java.nio.file.Path

import akka.actor.ActorSystem
import cromwell.backend.io.JobPaths
import cromwell.backend.{BackendJobDescriptorKey, BackendWorkflowDescriptor}
import cromwell.core.CallContext
import cromwell.services.metadata.CallMetadataKeys

object JesJobPaths {
  def apply(jobKey: BackendJobDescriptorKey, workflowDescriptor: BackendWorkflowDescriptor,
            jesConfiguration: JesConfiguration)(implicit actorSystem: ActorSystem): JesJobPaths = {
    new JesJobPaths(jobKey, workflowDescriptor, jesConfiguration)
  }

  val JesLogPathKey = "jesLog"
  val GcsExecPathKey = "gcsExec"
}

class JesJobPaths(val jobKey: BackendJobDescriptorKey, workflowDescriptor: BackendWorkflowDescriptor,
                  jesConfiguration: JesConfiguration)(implicit actorSystem: ActorSystem) extends
  JesWorkflowPaths(workflowDescriptor, jesConfiguration)(actorSystem) with JobPaths {

  val jesLogBasename = {
    val index = jobKey.index.map(s => s"-$s").getOrElse("")
    s"${jobKey.scope.unqualifiedName}$index"
  }

  override val returnCodeFilename: String = s"$jesLogBasename-rc.txt"
  override val stdoutFilename: String = s"$jesLogBasename-stdout.log"
  override val stderrFilename: String = s"$jesLogBasename-stderr.log"
  override val scriptFilename: String = "exec.sh"

  val jesLogFilename: String = s"$jesLogBasename.log"
  lazy val jesLogPath: Path = callExecutionRoot.resolve(jesLogFilename)
  
  lazy val callContext = CallContext(callExecutionRoot, stdoutFilename, stderrFilename)

  /*
  TODO: Move various monitoring files path generation here.

  "/cromwell_root" is a well known path, called in the regular JobPaths callDockerRoot.
  This JesCallPaths should know about that root, and be able to create the monitoring file paths.
  Instead of the AsyncActor creating the paths, the paths could then be shared with the CachingActor.

  Those monitoring paths could then be returned by metadataFiles and detritusFiles.
   */

  override lazy val customMetadataPaths = Map(
    CallMetadataKeys.BackendLogsPrefix + ":log" -> jesLogPath
  ) ++ (
    monitoringPath map { p => Map(JesMetadataKeys.MonitoringLog -> p) } getOrElse Map.empty  
  )

  override lazy val customDetritusPaths: Map[String, Path] = Map(
    JesJobPaths.GcsExecPathKey -> script,
    JesJobPaths.JesLogPathKey -> jesLogPath
  )
}
