package cromwell.backend.io

import com.typesafe.config.Config
import cromwell.backend.{BackendJobDescriptorKey, BackendWorkflowDescriptor}
import cromwell.core.path.PathFactory
import net.ceedubs.ficus.Ficus._
import cromwell.core.path.PathBuilder

class WorkflowPathsWithDocker(val workflowDescriptor: BackendWorkflowDescriptor,
                              val config: Config,
                              val pathBuilders: List[PathBuilder] = WorkflowPaths.DefaultPathBuilders) extends WorkflowPaths {

  val DockerRootString = config.as[Option[String]]("dockerRoot").getOrElse("/root")
  val DockerRoot = PathFactory.buildPath(DockerRootString, pathBuilders).toAbsolutePath
  val dockerWorkflowRoot = workflowPathBuilder(DockerRoot)

  override def toJobPaths(jobKey: BackendJobDescriptorKey,
                          jobWorkflowDescriptor: BackendWorkflowDescriptor): JobPathsWithDocker = {
    new JobPathsWithDocker(jobKey, jobWorkflowDescriptor, config, pathBuilders)
  }
}
