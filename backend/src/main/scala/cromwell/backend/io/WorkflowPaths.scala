package cromwell.backend.io

import java.nio.file.Path

import com.typesafe.config.Config
import cromwell.backend.{BackendJobDescriptorKey, BackendWorkflowDescriptor}
import cromwell.core.WorkflowOptions.FinalCallLogsDir
import cromwell.core.path.{DefaultPathBuilder, PathFactory}
import net.ceedubs.ficus.Ficus._

import scala.util.Try

object WorkflowPaths {
  val DefaultPathBuilders = List(DefaultPathBuilder)
}

trait WorkflowPaths extends PathFactory {
  def workflowDescriptor: BackendWorkflowDescriptor
  def config: Config
  
  protected lazy val executionRootString: String = config.as[Option[String]]("root").getOrElse("cromwell-executions")
  
  def getPath(url: String): Try[Path] = Try(PathFactory.buildPath(url, pathBuilders))
  
  // Rebuild potential intermediate call directories in case of a sub workflow
  protected def workflowPathBuilder(root: Path): Path = {
    workflowDescriptor.breadCrumbs.foldLeft(root)((acc, breadCrumb) => {
      breadCrumb.toPath(acc)
    }).resolve(workflowDescriptor.workflow.unqualifiedName).resolve(workflowDescriptor.id.toString + "/")
  }

  lazy val executionRoot: Path = PathFactory.buildPath(executionRootString, pathBuilders).toAbsolutePath
  lazy val workflowRoot: Path = workflowPathBuilder(executionRoot)
  lazy val finalCallLogsPath: Option[Path] =
    workflowDescriptor.getWorkflowOption(FinalCallLogsDir) map getPath map { _.get }

  def toJobPaths(jobKey: BackendJobDescriptorKey): JobPaths
}
