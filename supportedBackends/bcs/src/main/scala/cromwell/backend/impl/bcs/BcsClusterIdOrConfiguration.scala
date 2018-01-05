package cromwell.backend.impl.bcs


import scala.util.{Try, Success, Failure}
import scala.util.matching.Regex

final case class AutoClusterConfiguration(resourceType: String,
                                    instanceType: String,
                                    imageId: String,
                                    spotStrategy: Option[String] = None,
                                    spotPriceLimit: Option[Float] = None)


object BcsClusterIdOrConfiguration {
  type BcsClusterIdOrConfiguration = Either[String, AutoClusterConfiguration]

  val idPattern: Regex = """(cls-[^\s]+)""".r
  val resourceTypePattern = """(OnDemand|Spot)""".r
  val defaultResourceType = "OnDemand"

  val imageIdPattern = """([^\s]+)""".r

  val spotStrategyPattern = """(SpotWithPriceLimit|SpotAsPriceGo)""".r

  val spotPriceLimitPattern = """([01]\.\d{1,3})""".r

  // no suitable default instance type availabe
  val instanceTypePattern = """([be]cs[^\s]+)""".r

  val instanceAndImagePattern = s"""$instanceTypePattern\\s+$imageIdPattern""".r

  val resourceAndInstanceAndImagePattern = s"""$resourceTypePattern\\s+$instanceTypePattern\\s+$imageIdPattern""".r

  val spotPattern = s"""$resourceAndInstanceAndImagePattern\\s+$spotStrategyPattern\\s+$spotPriceLimitPattern""".r


  def parse(cluster: String): Try[BcsClusterIdOrConfiguration] = {
    cluster match {
      case idPattern(clusterId) => Success(Left(clusterId))
      case instanceAndImagePattern(instanceType, imageId) => Success(Right(AutoClusterConfiguration(defaultResourceType, instanceType, imageId)))
      case resourceAndInstanceAndImagePattern(resourceType, instanceType, imageId) => Success(Right(AutoClusterConfiguration(resourceType, instanceType, imageId)))
      case spotPattern(resourceType, instanceType, imageId, spotStrategy, spotPriceLimit) => Success(Right(AutoClusterConfiguration(resourceType, instanceType, imageId, Some(spotStrategy), Some(spotPriceLimit.toFloat))))
      case _ => Failure(new IllegalArgumentException("must be some string like 'cls-xxxx' or 'OnDemand ecs.s1.large img-ubuntu'"))
    }
  }
}
