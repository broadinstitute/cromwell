package centaur

import java.nio.file.Path

import org.scalatest.{ParallelTestExecution, Matchers, FlatSpec}

class CallCacheSpec extends FlatSpec with Matchers with ParallelTestExecution {

  def testCases(basePath: Path): List[WorkflowRequest] = {
    basePath.toFile.listFiles.toList collect { case x if x.isDirectory => x.toPath } map WorkflowRequest.apply
  }

  testCases(CentaurConfig.callCacheTestCasePath) foreach { case w =>
    if (w.name == "cacheWithinWF") {
      w.name should s"succesfully run ${w.name}" in {
        CacheFormulas.runCachingWorkflow(w).run.get
      }
    }
    if (w.name == "readFromCache" || w.name == "writeToCache" ) {
      w.name should s"successfully run ${w.name}" in {
        CacheFormulas.runCachingTurnedOffWorkflow(w).run.get
      }
    }
  }

  testCases(CentaurConfig.callCacheTestCasePath) foreach { case w =>
    if (w.name == "cacheWithinWF") {
      testCases(CentaurConfig.callCacheTestCasePath) foreach { case w2 =>
        if (w2.name == "cacheBetweenWF")
          w2.name should s"successfully run ${w.name}" in {
            CacheFormulas.runSequentialCachingWorkflow(w, w2).run.get
          }
      }
    }
  }

}