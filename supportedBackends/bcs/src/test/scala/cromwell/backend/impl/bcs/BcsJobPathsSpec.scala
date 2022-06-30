package cromwell.backend.impl.bcs

import common.mock.MockSugar
import cromwell.filesystems.oss.OssPath
import org.mockito.Mockito._

class BcsJobPathsSpec extends BcsTestUtilSpec with MockSugar {
  behavior of s"BcsJobPathsSpec"

  var root: OssPath = mockPathBuilder.build("oss://bcs-test/root/").getOrElse(throw new IllegalArgumentException())

  var workflowPath: BcsWorkflowPaths = {
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
