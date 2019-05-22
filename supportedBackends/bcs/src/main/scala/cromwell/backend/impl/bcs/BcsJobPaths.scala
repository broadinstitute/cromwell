package cromwell.backend.impl.bcs

import cromwell.backend.BackendJobDescriptorKey
import cromwell.backend.io.JobPaths
import cromwell.core.path.{DefaultPathBuilder, Path}

object BcsJobPaths {
	val BcsLogPathKey = "bcsLog"
	val BcsEnvExecKey = "exec"
	val BcsEnvCwdKey = "cwd"
	val BcsEnvStdoutKey = "stdout"
	val BcsEnvStderrKey = "stderr"
	val BcsCommandDirectory: Path = DefaultPathBuilder.get("/cromwell_root")
	val BcsTempInputDirectory: Path = DefaultPathBuilder.get("/cromwell_inputs")
	val BcsStdoutRedirectPath = "bcs-stdout"
	val BcsStderrRedirectPath = "bcs-stderr"
}

final case class BcsJobPaths(workflowPaths: BcsWorkflowPaths, jobKey: BackendJobDescriptorKey) extends JobPaths {

	import BcsJobPaths._

	val workerFileName = "worker"
	val workerPath = callRoot.resolve(workerFileName)
	val bcsStdoutPath = callRoot.resolve(BcsStdoutRedirectPath)
	val bcsStderrPath = callRoot.resolve(BcsStderrRedirectPath)
}
