package cromwell.binding

import cromwell.parser.BackendType
import cromwell.util.SampleWdl
import org.scalatest.prop.TableDrivenPropertyChecks._
import org.scalatest.{FlatSpec, Matchers}

class FqnResolverSpec extends FlatSpec with Matchers {
  val namespace = NamespaceWithWorkflow.load(SampleWdl.ScatterWdl.wdlSource(), BackendType.LOCAL)
  val callFqns = Table(
    "fqn",
    "w.A",
    "w.B",
    "w.C",
    "w.D",
    "w.E",
    "w.F",
    "w.G",
    "w.H"
  )

  val scatterFqns = Table(
    "fqn",
    "w.$scatter_0",
    "w.$scatter_0.$scatter_1",
    "w.$scatter_0.$scatter_2",
    "w.$scatter_3"
  )

  it should "be able to resolve a Workflow FQN" in {
    namespace.resolve("w") match {
      case Some(wf: Workflow) => wf.fullyQualifiedName shouldEqual "w"
      case Some(scope: Scope) => fail(s"FQN 'w' was expected to be a Workflow reference, but instead it references: $scope")
      case None => fail(s"Expecting to resolve FQN 'w' into a Workflow")
    }
  }

  forAll(callFqns) { (fqn: String) =>
    it should s"be able to resolve Call FQN: $fqn" in {
      namespace.resolve(fqn) match {
        case Some(call: Call) => call.fullyQualifiedName shouldEqual fqn
        case Some(scope: Scope) => fail(s"FQN '$fqn' was expected to be a Call reference, but instead it references: $scope")
        case None => fail(s"Expecting to resolve FQN '$fqn' into a Call")
      }
    }
  }

  forAll(scatterFqns) { (fqn: String) =>
    it should s"be able to resolve Scatter FQN: $fqn" in {
      namespace.resolve(fqn) match {
        case Some(scatter: Scatter) => scatter.fullyQualifiedNameWithIndexScopes shouldEqual fqn
        case Some(scope: Scope) => fail(s"FQN '$fqn' was expected to be a Scatter reference, but instead it references: $scope")
        case None => fail(s"Expecting to resolve FQN '$fqn' into a Scatter block")
      }
    }
  }
}
