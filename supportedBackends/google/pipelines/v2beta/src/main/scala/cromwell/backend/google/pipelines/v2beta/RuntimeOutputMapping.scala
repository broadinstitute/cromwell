package cromwell.backend.google.pipelines.v2beta

import common.util.StringUtil._
import cromwell.backend.google.pipelines.common.io.PipelinesApiWorkingDisk
import cromwell.core.path.Path

object RuntimeOutputMapping {

  /**
    * List of prefixes to be stripped away from runtime output paths before
    * appending them to the cloud call root to generate the delocalization path.
    * The goal is to reduce unnecessary long paths which keep repeating the workflow root
    * as the workflow progresses
    * 
    * For instance:
    * 
    * file:///cromwell_root/bucket/workflow_name/6d777414-5ee7-4c60-8b9e-a02ec44c398e/call-A/file.txt
    *  ->
    * call-A/file.txt
    * 
    * Which will be delocalized to
    * gs://bucket/workflow_name/6d777414-5ee7-4c60-8b9e-a02ec44c398e/call-B/call-A/file.txt
    *  instead of
    * gs://bucket/workflow_name/6d777414-5ee7-4c60-8b9e-a02ec44c398e/call-B/cromwell_root/bucket/workflow_name/6d777414-5ee7-4c60-8b9e-a02ec44c398e/call-A/file.txt
    */
  def prefixFilters(workflowRoot: Path): List[String] = List(
    "file://",
    "/",
    PipelinesApiWorkingDisk.MountPoint.pathAsString.relativeDirectory,
    workflowRoot.pathWithoutScheme.relativeDirectory
  )
}
