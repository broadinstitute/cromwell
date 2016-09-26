package centaur

import java.nio.file.Paths

import cats.Apply
import cats.data.Validated.{Invalid, Valid}
import cats.std.list._
import centaur.test.ErrorOr
import centaur.test.formulas.TestFormulas
import centaur.test.workflow.Workflow
import org.scalatest.{FlatSpec, Matchers, ParallelTestExecution}

object CallCacheSpec {
  val CallCachingWorkflowDir = Paths.get("src/main/resources/callCachingWorkflows")
  val ReadFromCacheTest = CallCachingWorkflowDir.resolve("readFromCache.test")
  val WriteToCacheTest = CallCachingWorkflowDir.resolve("writeToCache.test")
  val CacheWithinWf = CallCachingWorkflowDir.resolve("cacheWithinWF.test")
  val CacheBetweenWf = CallCachingWorkflowDir.resolve("cacheBetweenWf.test")
}

class CallCacheSpec extends FlatSpec with Matchers with ParallelTestExecution {
  import CallCacheSpec._

  "readFromCacheOff" should "not use call cache reading" in {
    Workflow.fromPath(ReadFromCacheTest) match {
      case Valid(w) => TestFormulas.runCachingTurnedOffWorkflow(w).run.get
      case Invalid(e) => fail(s"Could not read readFromCache test:\n -${e.unwrap.mkString("\n-")}")
    }
  }

  "cacheBetweenWf" should "successfully call cache between two workflows" in {
    Workflow.fromPath(CacheBetweenWf) match {
      case Valid(w) => TestFormulas.runSequentialCachingWorkflows(w, w)
      case Invalid(e) => fail(s"Could not read cacheWithinWf test:\n - ${e.unwrap.mkString("\n- ")}")
    }
  }
}
