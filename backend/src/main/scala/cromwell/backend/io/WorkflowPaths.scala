package cromwell.backend.io

import com.typesafe.config.Config
import cromwell.backend.{BackendJobDescriptor, BackendJobDescriptorKey, BackendWorkflowDescriptor}
import cromwell.core.WorkflowOptions.FinalCallLogsDir
import cromwell.core.path.{DefaultPathBuilder, Path, PathFactory}
import net.ceedubs.ficus.Ficus._

import scala.util.Try

object WorkflowPaths {
  val DefaultPathBuilders = List(DefaultPathBuilder)
  val DefaultExecutionRootString: String = "cromwell-executions"
}

trait WorkflowPaths extends PathFactory {
  def workflowDescriptor: BackendWorkflowDescriptor
  def config: Config

  /**
    * Path (as a String) of the root directory Cromwell should use for ALL workflows.
    */
  protected lazy val executionRootString: String =
    config.as[Option[String]]("root").getOrElse(WorkflowPaths.DefaultExecutionRootString)

  /**
    * Implementers of this trait might override this to provide an appropriate prefix corresponding to the execution root
    * of the current workflow. This prefix would be prepended to the list of prefixes provided in workflow options for
    * searching cache hits.
    */
  lazy val callCacheRootPrefix: Option[String] = None

  /**
    * Path of the root directory Cromwell will use for ALL workflows.
    */
  lazy val executionRoot: Path = buildPath(executionRootString).toAbsolutePath

  /**
    * This MUST be a directory that can be accessed by Cromwell in a read and write fashion.
    * It will contain all files that Cromwell will create for this workflow as well as outputs
    * from the jobs.
    * The structure under this directory is relatively stable but should not be relied upon too heavily
    * as it might be subject to changes.
    *
    * This can be either a Shared Filesystem or Cloud Filesystem Path as long as Cromwell's filesystems are configured accordingly.
    */
  lazy val workflowRoot: Path = workflowPathBuilder(executionRoot)

  /**
    * Attempts to build a cromwell.core.Path from the String using the available Filesystems.
    */
  def getPath(url: String): Try[Path] = Try(buildPath(url))

  // Rebuild potential intermediate call directories in case of a sub workflow
  protected def workflowPathBuilder(root: Path): Path =
    workflowDescriptor.breadCrumbs
      .foldLeft(root)((acc, breadCrumb) => breadCrumb.toPath(acc))
      .resolve(workflowDescriptor.callable.name)
      .resolve(workflowDescriptor.id.toString + "/")

  lazy val finalCallLogsPath: Option[Path] =
    workflowDescriptor.getWorkflowOption(FinalCallLogsDir) map getPath map { _.get }

  def toJobPaths(jobDescriptor: BackendJobDescriptor): JobPaths =
    toJobPaths(jobDescriptor.key, jobDescriptor.workflowDescriptor)

  /**
    * Creates job paths using the key and workflow descriptor.
    *
    * NOTE: For sub workflows, the jobWorkflowDescriptor will be different than the WorkflowPaths.workflowDescriptor.
    *
    * @param jobKey                The key for the job.
    * @param jobWorkflowDescriptor The workflow descriptor for the job.
    * @return The paths for the job.
    */
  def toJobPaths(jobKey: BackendJobDescriptorKey, jobWorkflowDescriptor: BackendWorkflowDescriptor): JobPaths =
    // If the descriptors are the same, no need to create a new WorkflowPaths
    if (workflowDescriptor == jobWorkflowDescriptor) toJobPaths(this, jobKey)
    else toJobPaths(withDescriptor(jobWorkflowDescriptor), jobKey)

  protected def toJobPaths(workflowPaths: WorkflowPaths, jobKey: BackendJobDescriptorKey): JobPaths

  protected def withDescriptor(workflowDescriptor: BackendWorkflowDescriptor): WorkflowPaths
}
