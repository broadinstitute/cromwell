package cromwell.engine.backend

package object pbs {
  /**
    * PBE This class is a temporary workaround for the fact that Backends are not call-scoped yet.
    * When they are, all those values can be declared directly in the backend, as it will have a BackendCallJobDescriptor parameter.
    * All values have been turned into defs because Value Classes only support defs.
    */
  implicit class PbsJobDescriptor(val jobDescriptor: BackendCallJobDescriptor) extends AnyVal {
    def returnCode = jobDescriptor.callRootPath.resolve("rc")
    def stdout = jobDescriptor.callRootPath.resolve("stdout")
    def stderr = jobDescriptor.callRootPath.resolve("stderr")
    def script = jobDescriptor.callRootPath.resolve("script.sh")
  }
}
