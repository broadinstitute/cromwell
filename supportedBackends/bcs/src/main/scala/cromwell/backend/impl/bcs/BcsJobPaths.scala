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
}

final case class BcsJobPaths(workflowPaths: BcsWorkflowPaths, jobKey: BackendJobDescriptorKey) extends JobPaths {

	// alibaba cloud's batchcompute service can only support tar.gz formatted package.
	val workerFileName = "worker.tar.gz"
	val workerPath = workflowPaths.workflowRoot.resolve(workerFileName)
}
