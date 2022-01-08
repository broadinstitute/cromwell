package cromwell.backend.impl.bcs

import org.scalatest.prop.Tables.Table
import org.scalatest.prop.TableDrivenPropertyChecks._
import org.scalatest.TryValues._

import scala.util.Failure


class BcsClusterIdOrConfigurationSpec extends BcsTestUtilSpec {
  behavior of s"BcsClusterIdOrConfiguration"

  val clusterIdTable = Table(
    ("unparsed", "parsed"),
    ("cls-xxxx", Option("cls-xxxx")),
    ("job-xxxx", None)
  )

  it should "parse correct cluster id" in {
    forAll(clusterIdTable) { (unparsed, parsed) =>
      BcsClusterIdOrConfiguration.idPattern.findFirstIn(unparsed) shouldEqual(parsed)
    }
  }

  val resourceTypeTable = Table(
    ("unparsed", "parsed"),
    ("OnDemand", Option("OnDemand")),
    ("Spot", Option("Spot")),
    ("Other", None)
  )

  it should "parse correct resource type" in {
    forAll(resourceTypeTable) { (unparsed, parsed) =>
      BcsClusterIdOrConfiguration.resourceTypePattern.findFirstIn(unparsed) shouldEqual(parsed)
    }
  }


  val instanceTypeTable = Table(
    ("unparsed", "parsed"),
    ("ecs.s1.large", Option("ecs.s1.large")),
    ("bcs.s1.large", Option("bcs.s1.large"))
  )

  it should "parse correct instance type" in {
    forAll(instanceTypeTable) { (unparsed, parsed) =>
      BcsClusterIdOrConfiguration.instanceTypePattern.findFirstIn(unparsed) shouldEqual parsed
    }
  }

  val spotStrategyTable = Table(
    ("unparsed", "parsed"),
    ("SpotWithPriceLimit", Option("SpotWithPriceLimit")),
    ("SpotAsPriceGo", Option("SpotAsPriceGo"))
  )


  it should "parse correct spot strategy" in {
    forAll(spotStrategyTable) { (unparsed, parsed) =>
      BcsClusterIdOrConfiguration.spotStrategyPattern.findFirstIn(unparsed) shouldEqual parsed
    }
  }

  val spotPriceLimitTable = Table(
    ("unparsed", "parsed"),
    ("1.0", Option(1.0.toFloat)),
    ("0.1", Option(0.1.toFloat)),
    ("0.12", Option(0.12.toFloat)),
    ("0.123", Option(0.123.toFloat))
  )

  it should "parse correct spot price limit" in {
    forAll(spotPriceLimitTable) { (unparsed, parsed) =>
      BcsClusterIdOrConfiguration.spotPriceLimitPattern.findFirstIn(unparsed) map {limit => limit.toFloat} shouldEqual parsed
    }
  }

  val validClusterInfoTable = Table(
    ("unparsed", "parsed"),
    ("cls-id", Left("cls-id")),
    ("OnDemand ecs.s1.large img-test", Right(AutoClusterConfiguration("OnDemand", "ecs.s1.large", "img-test"))),
    ("OnDemand ecs.s1.large img-test cls-test", Right(AutoClusterConfiguration("OnDemand", "ecs.s1.large", "img-test", clusterId = Option("cls-test")))),
    ("ecs.s1.large img-test", Right(AutoClusterConfiguration("OnDemand", "ecs.s1.large", "img-test"))),
    ("ecs.s1.large img-test cls-test", Right(AutoClusterConfiguration("OnDemand", "ecs.s1.large", "img-test", clusterId = Option("cls-test")))),
    ("Spot ecs.s1.large img-test SpotWithPriceLimit 0.1", Right(AutoClusterConfiguration("Spot", "ecs.s1.large", "img-test", Option("SpotWithPriceLimit"), Option(0.1.toFloat)))),
    ("Spot ecs.s1.large img-test SpotWithPriceLimit 0.1 cls-test", Right(AutoClusterConfiguration("Spot", "ecs.s1.large", "img-test", Option("SpotWithPriceLimit"), Option(0.1.toFloat), Option("cls-test")))),
    ("Spot ecs.s1.large img-test SpotAsPriceGo 0.1", Right(AutoClusterConfiguration("Spot", "ecs.s1.large", "img-test", Option("SpotAsPriceGo"), Option(0.1.toFloat)))),
    ("Spot ecs.s1.large img-test SpotAsPriceGo 0.1 cls-test", Right(AutoClusterConfiguration("Spot", "ecs.s1.large", "img-test", Option("SpotAsPriceGo"), Option(0.1.toFloat), Option("cls-test")))),

  )


  it should "parse correct cluster id or cluster configuration" in {
    forAll(validClusterInfoTable) { (unparsed, parsed) =>
      BcsClusterIdOrConfiguration.parse(unparsed).success.value shouldEqual parsed
    }
  }

  val invalidClusterInfos = List("OnDemand", "", "OnDemand other", "other ecs.s1.large")

  invalidClusterInfos foreach { unparsed =>
    it should s"throw when parsing $unparsed" in {
      BcsClusterIdOrConfiguration.parse(unparsed) shouldBe a [Failure[_]]
    }
  }
}
