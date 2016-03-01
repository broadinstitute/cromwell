package centaur

import java.nio.file.Path

import scala.language.postfixOps
import org.scalatest.{Matchers, FlatSpec, ParallelTestExecution}

class CentaurSpec extends FlatSpec with Matchers with ParallelTestExecution {
  /**
    * FIXME: There's a similar comment in application.conf. When we start supporting behavior beyond "run everything
    * and make sure it works" either there'll need to be multiple test case directories or some other metadata
    * to provide classification on how the workflow is supposed to be handled
    */
  def testCases(basePath: Path): List[WorkflowRequest] = {
    basePath.toFile.listFiles.toList collect { case x if x.isDirectory => x.toPath } map WorkflowRequest.apply
  }


  testCases(CentaurConfig.successfulTestCasePath) foreach { case w =>
    w.name should "successfully run" in {
      TestFormulas.runSuccessfulWorkflow(w).run.get
      Thread.sleep(1000)
    }
  }

  testCases(CentaurConfig.failingTestCasePath) foreach { case w =>
    w.name should "fail" in {
      TestFormulas.runFailingWorkflow(w).run.get
      Thread.sleep(1000)
    }
  }
}
