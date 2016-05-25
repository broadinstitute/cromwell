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
  val CacheWithinWf = CallCachingWorkflowDir.resolve("cacheWithinWf.test")
  val CacheBetweenWf = CallCachingWorkflowDir.resolve("cacheBetweenWF.test")
}

class CallCacheSpec extends FlatSpec with Matchers with ParallelTestExecution {
  import CallCacheSpec._

  "readFromCache" should "not use call caching" in {
    Workflow.fromPath(ReadFromCacheTest) match {
      case Valid(w) => TestFormulas.runCachingTurnedOffWorkflow(w).run.get
      case Invalid(e) => fail(s"Could not read readFromCache test:\n -${e.unwrap.mkString("\n-")}")
    }
  }

  "writeToCache" should "not use call caching" in {
    Workflow.fromPath(WriteToCacheTest) match {
      case Valid(w) => TestFormulas.runCachingTurnedOffWorkflow(w).run.get
      case Invalid(e) => fail(s"Could not read writeToCache test:\n - ${e.unwrap.mkString("\n- ")}")
    }
  }

  "cacheWithinWf" should "use call caching" in {
    Workflow.fromPath(CacheWithinWf) match {
      case Valid(w) => TestFormulas.runSuccessfulWorkflowAndVerifyMetadata(w).run.get
      case Invalid(e) => fail(s"Could not read cacheWithinWf test:\n - ${e.unwrap.mkString("\n- ")}")
    }
  }

  "cacheBetweenWf" should "use call caching" in {
    val cacheWithinWf = Workflow.fromPath(CacheWithinWf)
    val cacheBetweenWf = Workflow.fromPath(CacheBetweenWf)
    Apply[ErrorOr].map2(cacheWithinWf, cacheBetweenWf)((w, b) => TestFormulas.runSequentialCachingWorkflow(w, b)) match {
      case Valid(t) => t.run.get
      case Invalid(e) => fail(s"Could not build workflows for cacheBetwenWf test:\n - ${e.unwrap.mkString("\n- ")}")
    }
  }
}
