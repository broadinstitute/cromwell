package centaur

import java.nio.file.Path

import scala.language.postfixOps
import org.scalatest.{Tag, Matchers, FlatSpec, ParallelTestExecution}

/**
  * Tests with this tag expect their initial submission to Cromwell to fail
  */
object SubmissionFailureTest extends Tag("SubmissionFailureTest")
/**
  * Tests with this tag expect to fail
  */
object FailureTest extends Tag("FailureTest")

class CentaurSpec extends FlatSpec with Matchers with ParallelTestExecution {

  def testCases(basePath: Path): List[WorkflowRequest] = {
    basePath.toFile.listFiles.toList collect { case x if x.isDirectory => x.toPath } map WorkflowRequest.apply
  }

  testCases(CentaurConfig.successfulTestCasePath) foreach { case w =>
    w.name should "successfully run ${w.name}" in {
      TestFormulas.runSuccessfulWorkflow(w).run.get
      Thread.sleep(1000)
    }
  }

  testCases(CentaurConfig.failingTestCasePath) foreach { case w =>
    w.name should s"fail ${w.name}" taggedAs FailureTest in {
      TestFormulas.runFailingWorkflow(w).run.get
      Thread.sleep(1000)
    }
  }

  testCases(CentaurConfig.submissionFailureTestCasePath) foreach { case w =>
    w.name should s"fail ${w.name}" taggedAs SubmissionFailureTest in {
      // TODO: This returns a string error message. With extra test metadata, we can verify that the error message is correct!
      TestFormulas.runSubmissionFailureWorkflow(w).run.get
      Thread.sleep(1000)
    }
  }
}
