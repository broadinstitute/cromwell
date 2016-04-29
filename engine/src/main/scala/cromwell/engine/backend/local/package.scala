package cromwell.engine.backend

package object local {
  /**
    * PBE This class is a temporary workaround for the fact that Backends are not call-scoped yet.
    * When they are, all those values can be declared directly in the backend, as it will have a BackendCallJobDescriptor parameter.
    * All values have been turned into defs because Value Classes only support defs.
    */
  implicit class LocalJobDescriptor(val jobDescriptor: OldStyleBackendCallJobDescriptor) extends AnyVal {
    def containerCallRoot = jobDescriptor.callRuntimeAttributes.docker match {
      case Some(docker) => jobDescriptor.callRootPathWithBaseRoot(OldStyleLocalBackend.ContainerRoot)
      case None => jobDescriptor.callRootPath
    }
    def returnCode = jobDescriptor.callRootPath.resolve("rc")
    def stdout = jobDescriptor.callRootPath.resolve("stdout")
    def stderr = jobDescriptor.callRootPath.resolve("stderr")
    def script = jobDescriptor.callRootPath.resolve("script")
  }
}
