package cromwell.webservice

import cromwell.CromwellTestkitSpec
import cromwell.engine.db.ExecutionDatabaseKey

import scalaz._

object CallCachingParametersSpec {
  val CallNames = Seq("three_step.cgrep", "scatter.B")
  val AllowTrue = QueryParameter("allow", "true")
  val AllowFalse = QueryParameter("allow", "false")
}

class CallCachingParametersSpec extends CromwellTestkitSpec("CallCachingParametersSpec") {

  import CallCachingParametersSpec._

  "CallCachingParameters call name validation" should {
    "be accepted if empty" in {
      CallCachingParameters.validateCallName(None) match {
        case Success(k) => k shouldBe None
        case Failure(x) => fail("Unexpected failure: " + x.list.mkString(", "))
      }
    }

    "be accepted if not indexed" in {
      CallNames foreach { name =>
        CallCachingParameters.validateCallName(Option(name)) match {
          case Success(k) => k shouldEqual Some(ExecutionDatabaseKey(name, None))
          case Failure(x) => fail("Unexpected failure: " + x.list.mkString(", "))
        }
      }
    }

    "be accepted if indexed" in {
      CallNames foreach { name =>
        CallCachingParameters.validateCallName(Option(name + ".0")) match {
          case Success(k) => k shouldEqual Some(ExecutionDatabaseKey(name, Some(0)))
          case Failure(x) => fail("Unexpected failure: " + x.list.mkString(", "))
        }
      }
    }

    "be rejected if a scatter" in {
      CallCachingParameters.validateCallName(Option("scattergather.$scatter_0")) match {
        case Success(k) => fail("Unexpected Success: " + k)
        case Failure(_) =>
      }
    }
  }

  "CallCachingParameters recognized keys validation" should {
    "disallow unrecognized keys" in {
      CallCachingParameters.validateRecognizedKeys(Seq(AllowTrue, QueryParameter("who", "dat"))) match {
        case Success(k) => fail("Unexpected Success: " + k)
        case Failure(_) =>
      }
    }
  }

  "CallCachingParameters allow validation" should {
    "disallow non-boolean values" in {
      CallCachingParameters.validateAllow(Seq(QueryParameter("allow", "sure why not"))) match {
        case Success(k) => fail("Unexpected Success: " + k)
        case Failure(_) =>
      }
    }

    // Arguably this shouldn't allow multiple `allow`s at all, but for now the validation lets them slide
    // as long as they're all the same.
    "disallow incoherent settings" in {
      CallCachingParameters.validateAllow(Seq(AllowTrue, AllowFalse)) match {
        case Success(k) => fail("Unexpected Success: " + k)
        case Failure(_) =>
      }
    }
  }
}
