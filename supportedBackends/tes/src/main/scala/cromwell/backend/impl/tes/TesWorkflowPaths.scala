package cromwell.backend.impl.tes

import com.typesafe.config.Config
import cromwell.backend.io.WorkflowPaths
import cromwell.backend.{BackendJobDescriptorKey, BackendWorkflowDescriptor}
import cromwell.core.path.{PathBuilder, PathFactory}
import net.ceedubs.ficus.Ficus._

case class TesWorkflowPaths(override val workflowDescriptor: BackendWorkflowDescriptor,
                       override val config: Config,
                       override val pathBuilders: List[PathBuilder] = WorkflowPaths.DefaultPathBuilders) extends WorkflowPaths {

  val DockerRootString = config.as[Option[String]]("dockerRoot").getOrElse("/cromwell-executions")
  var DockerRoot = PathFactory.buildPath(DockerRootString, pathBuilders)
  if (!DockerRoot.isAbsolute) {
    DockerRoot = PathFactory.buildPath("/".concat(DockerRootString), pathBuilders)
  }
  val dockerWorkflowRoot = workflowPathBuilder(DockerRoot)

  override def toJobPaths(workflowPaths: WorkflowPaths,
                          jobKey: BackendJobDescriptorKey): TesJobPaths = {
    new TesJobPaths(workflowPaths.asInstanceOf[TesWorkflowPaths], jobKey)
  }

  override protected def withDescriptor(workflowDescriptor: BackendWorkflowDescriptor): WorkflowPaths = this.copy(workflowDescriptor = workflowDescriptor)
}
