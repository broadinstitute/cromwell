package cromwell.backend.impl.bcs

import cromwell.backend.BackendJobDescriptorKey
import cromwell.filesystems.oss.OssPath
import org.mockito.Mockito.when
import wdl4s.wdl.WdlTaskCall

object BcsJobPathsSpec {
}


class BcsJobPathsSpec extends BcsTestUtilSpec {
  behavior of s"BcsJobPathsSpec"

  var root: OssPath = mockPathBuiler.build("oss://bcs-test/root/").getOrElse(throw new IllegalArgumentException())

  var workflowPath = {
    val workflowPaths = mock[BcsWorkflowPaths]

    when(workflowPaths.workflowRoot).thenReturn(root)
    workflowPaths
  }

  def name = "test"

  override val jobKey = {
    val key = mock[BackendJobDescriptorKey]
    val scope = mock[WdlTaskCall]
    when(scope.unqualifiedName).thenReturn(name)
    when(key.attempt).thenReturn(0)
    when(key.index).thenReturn(None)

    when(key.scope).thenReturn(scope)

    key
  }

  it should "have right package name" in {
    val jobPath = BcsJobPaths(workflowPath, jobKey)
    jobPath.workerPath shouldEqual root.resolve(s"${jobPath.workerFileName}")
  }

}
