package cromwell.backend.impl.bcs

import cromwell.filesystems.oss.OssPath
import org.mockito.Mockito.when

class BcsJobPathsSpec extends BcsTestUtilSpec {
  behavior of s"BcsJobPathsSpec"

  var root: OssPath = mockPathBuilder.build("oss://bcs-test/root/").getOrElse(throw new IllegalArgumentException())

  var workflowPath = {
    val workflowPaths = mock[BcsWorkflowPaths]

    when(workflowPaths.workflowRoot).thenReturn(root)
    workflowPaths
  }

  def name = "test"

  it should "have right package name" in {
    val jobPath = BcsJobPaths(workflowPath, jobKey)
    jobPath.workerPath shouldEqual jobPath.callRoot.resolve(s"${jobPath.workerFileName}")
  }
}
