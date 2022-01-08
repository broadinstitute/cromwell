package cromwell.backend.impl.bcs

import com.typesafe.config.ConfigValueFactory
import cromwell.backend.impl.bcs.callcaching.UseOriginalCachedOutputs

class BcsConfigurationSpec extends BcsTestUtilSpec {
  behavior of "BcsConfiguration"
  type ValueOrDelete = Either[Boolean, AnyRef]

  def backendConfiguration = BcsTestUtilSpec.BcsBackendConfigurationDescriptor
  def defaultBackendConfig = BcsTestUtilSpec.BcsBackendConfigurationDescriptor.backendConfig

  it should "have correct bcs region" in {
    val region = "cn-hangzhou"
    val configs = Map("region" -> Right(region))
    val conf = withConfig(configs)
    conf.bcsRegion shouldEqual Some(region)
  }

  it should "have correct bcs access id and key" in {
    val id = "test-access-id"
    val key = "test-access-key"
    val configs = Map("access-id" -> Right(id), "access-key" -> Right(key))
    val conf = withConfig(configs)
    conf.bcsAccessId shouldEqual Some(id)
    conf.bcsAccessKey shouldEqual Some(key)
  }

  it should "have correct bcs callcaching strategy" in {
    val region = "cn-hangzhou"
    val configs = Map("region" -> Right(region))
    val conf = withConfig(configs)
    conf.duplicationStrategy shouldEqual UseOriginalCachedOutputs
  }


  private def withConfig(configs: Map[String, ValueOrDelete]) = {
    var descriptor = BcsTestUtilSpec.BcsBackendConfigurationDescriptor.copy()
    for ((key, value) <- configs) {
      value match {
        case Left(_) => descriptor = BcsTestUtilSpec.BcsBackendConfigurationDescriptor.copy(backendConfig = descriptor.backendConfig.withoutPath(key))
        case Right(v) => descriptor = BcsTestUtilSpec.BcsBackendConfigurationDescriptor.copy(backendConfig = descriptor.backendConfig.withValue(key, ConfigValueFactory.fromAnyRef(v)))
      }
    }
    new BcsConfiguration(descriptor)
  }

}
