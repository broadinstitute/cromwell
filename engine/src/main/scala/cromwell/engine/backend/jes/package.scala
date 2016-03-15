package cromwell.engine.backend

import java.nio.file.Paths

import wdl4s.values.WdlFile

import scala.util.Try

package object jes {
  /**
    * PBE This class is a temporary workaround for the fact that Backends are not call-scoped yet.
    * When they are, all those values can be declared directly in the backend, as it will have a BackendCallJobDescriptor parameter.
    * All values have been turned into defs because Value Classes only support defs.
    */
  implicit class JesJobDescriptor(val jobDescriptor: BackendCallJobDescriptor) extends AnyVal {
    import JesBackend._
    import better.files._

    // TODO: Assuming that runtimeAttributes.disks always has a 'local-disk'
    def workingDisk = jobDescriptor.callRuntimeAttributes.disks.find(_.name == JesWorkingDisk.Name).get
    def gcsExecPath = jobDescriptor.callRootPath.resolve(JesExecScript)
    def defaultMonitoringOutputPath = jobDescriptor.callRootPath.resolve(JesMonitoringLogFile)
    def returnCodeFilename = jesReturnCodeFilename(jobDescriptor.key)
    def jesStdoutGcsPath = jobDescriptor.callRootPath.resolve(jesLogStdoutFilename(jobDescriptor.key))
    def jesStderrGcsPath = jobDescriptor.callRootPath.resolve(jesLogStderrFilename(jobDescriptor.key))
    def jesLogGcsPath = jobDescriptor.callRootPath.resolve(jesLogFilename(jobDescriptor.key))
    def returnCodeGcsPath = jobDescriptor.callRootPath.resolve(returnCodeFilename)
    def rcJesOutput = JesFileOutput(returnCodeFilename, returnCodeGcsPath.toString, Paths.get(returnCodeFilename), workingDisk)
    def cmdInput = JesFileInput(ExecParamName, gcsExecPath.toString, Paths.get(JesExecScript), workingDisk)
    def maxPreemption = jobDescriptor.callRuntimeAttributes.preemptible
    def jesCommandLine = s"/bin/bash ${cmdInput.containerPath.toAbsolutePath.toString}"
    def standardParameters = Seq(rcJesOutput)
    def stdoutStderr: CallLogs = {
      CallLogs(
        stdout = WdlFile(jesStdoutGcsPath.toString),
        stderr = WdlFile(jesStderrGcsPath.toString),
        Option(Map("log" -> WdlFile(jesLogGcsPath.toString)))
      )
    }
    /**
      * Determines the maximum number of times a call can be started with a Preemptible VM.
      * TODO: Use configuration as a way to set this globally.
      * Currently workflow options act as default for runtime attributes, configuration could do the same for workflow options.
      */
    def preemptible = JesBackend.preemptible(jobDescriptor, maxPreemption)

    def downloadRcFile = Try(returnCodeGcsPath.toAbsolutePath.contentAsString)
    /**
      * Determine the output directory for the files matching a particular glob.
      */
    def globOutputPath(glob: String) = jobDescriptor.callRootPath.resolve(globDirectory(glob)).toString
  }
}
