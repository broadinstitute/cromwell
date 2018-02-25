package cromwell.backend.impl.bcs

class BcsWorkflowPathsSpec extends BcsTestUtilSpec {

  behavior of s"BcsWorkflowPaths"

  it should "have correct input workflow mapping" in {

    import BcsTestUtilSpec._

    val paths = BcsWorkflowPaths(workflowDescriptor, BcsBackendConfig, mockPathBuilders)

    val workflowInput = paths.getWorkflowInputMounts
    workflowInput shouldBe a[BcsInputMount]
    workflowInput.src shouldEqual(paths.workflowRoot)
    workflowInput.dest.pathAsString.startsWith(BcsJobPaths.BcsTempInputDirectory.pathAsString) shouldBe true
    // DefaultPathBuilder always remove ending '/' from directory path.
    workflowInput.dest.pathAsString.endsWith(paths.workflowRoot.pathWithoutScheme.stripSuffix("/")) shouldBe true
  }

  it should "have correct job paths" in  {

    import BcsTestUtilSpec._

    val paths = BcsWorkflowPaths(workflowDescriptor, BcsBackendConfig, mockPathBuilders)

    val jobPaths = paths.toJobPaths(paths, jobKey)
    jobPaths shouldBe a [BcsJobPaths]
  }

}
