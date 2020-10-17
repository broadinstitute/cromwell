package cromwell.engine.workflow.lifecycle

import com.typesafe.config.ConfigFactory
import common.assertion.CromwellTimeoutSpec
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import org.scalatest.prop.TableDrivenPropertyChecks

import scala.util.{Failure, Success, Try}

class ValidatingCachingConfigSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers with TableDrivenPropertyChecks {

    it should "run config tests" in {
      val cases = Table(
        ("config"                 , "exceptionMessage"                                                              ),
        ("enabled = not-a-boolean", "String: 1: enabled has type STRING rather than BOOLEAN"                        ),
        ("enabled = true"         , true                                                                            ),
        ("enabled = false"        , false                                                                           ),
        ("enabled = 1"            , "String: 1: enabled has type NUMBER rather than BOOLEAN"                        ),
        (""                       , "String: 1: No configuration setting found for key 'enabled'"                   )
      )

      forEvery(cases) { (config, expected) =>
        val rootConfig = ConfigFactory.parseString(config)
        Try(rootConfig.getBoolean("enabled")) match {
          case Success(what) => what shouldBe a [java.lang.Boolean]
          case Failure(exception) => exception.getMessage should be (expected)
        }
      }
    }
  }


