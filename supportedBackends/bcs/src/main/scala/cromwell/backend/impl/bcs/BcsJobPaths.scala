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

final case class BcsJobPaths(override val workflowPaths: BcsWorkflowPaths, jobKey: BackendJobDescriptorKey) extends JobPaths {

	val baseName = {
		val index = jobKey.index.map(s => s"-$s").getOrElse("")
		s"${jobKey.scope.unqualifiedName}$index"
	}

	// alibaba cloud's batchcompute service can only support tar.gz formatted package.
	val workerFileName = "worker.tar.gz"
	val worker = workflowPaths.workflowRoot.resolve(workerFileName)
}
