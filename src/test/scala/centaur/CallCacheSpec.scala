package centaur

import java.nio.file.Paths

import cats.data.Validated.{Invalid, Valid}
import centaur.test.formulas.TestFormulas
import centaur.test.workflow.Workflow
import org.scalatest.{FlatSpec, Matchers, ParallelTestExecution}

object CallCacheSpec {
  val CallCachingWorkflowDir = Paths.get("src/main/resources/callCachingWorkflows")
  val ReadFromCacheTest = CallCachingWorkflowDir.resolve("readFromCache.test")
  val FloatingTagsTest = CallCachingWorkflowDir.resolve("floatingTags.test")
  val WriteToCacheTest = CallCachingWorkflowDir.resolve("writeToCache.test")
  val CacheWithinWf = CallCachingWorkflowDir.resolve("cacheWithinWF.test")
  val CacheBetweenWf = CallCachingWorkflowDir.resolve("cacheBetweenWf.test")
}

class CallCacheSpec extends FlatSpec with Matchers with ParallelTestExecution {
  import CallCacheSpec._

  "readFromCacheOff" should "not use call cache reading" in {
    Workflow.fromPath(ReadFromCacheTest) match {
      case Valid(w) => TestFormulas.runCachingTurnedOffWorkflow(w).run.get
      case Invalid(e) => fail(s"Could not read readFromCache test:\n -${e.toList.mkString("\n-")}")
    }
  }

  "cacheBetweenWf" should "successfully call cache between two workflows" in {
    Workflow.fromPath(CacheBetweenWf) match {
      case Valid(w) => TestFormulas.runSequentialCachingWorkflows(w, w).run.get
      case Invalid(e) => fail(s"Could not read cacheWithinWf test:\n - ${e.toList.mkString("\n- ")}")
    }
  }

  "floatingTags" should "not hit the cache when tasks have floating docker tags" in {
    Workflow.fromPath(FloatingTagsTest) match {
      case Valid(w) => TestFormulas.runCachingTurnedOffWorkflow(w).run.get
      case Invalid(e) => fail(s"Could not read floatingTags test:\n -${e.toList.mkString("\n-")}")
    }
  }
}
