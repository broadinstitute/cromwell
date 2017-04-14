package centaur

import cats.data.Validated.Valid
import centaur.api.CromwellBackendsCompanion
import centaur.test.formulas.TestFormulas
import centaur.test.workflow.Workflow
import org.scalatest.FlatSpec

class BackendNoDockerSpec extends FlatSpec {

  val backendWithNoDockerVerb = "backendWithNoDocker" should "ignore docker attribute and call cache even with tag"
  val BackendWithNoDocker = CallCacheSpec.CallCachingWorkflowDir.resolve("backendWithNoDocker.test")
  
  Workflow.fromPath(BackendWithNoDocker) match {
    case Valid(w) if w.backends forall CromwellBackendsCompanion.supportedBackends.contains =>
      backendWithNoDockerVerb in { TestFormulas.runSequentialCachingWorkflows(w, w).run.get }
    case _ => backendWithNoDockerVerb ignore { }
  }
}
