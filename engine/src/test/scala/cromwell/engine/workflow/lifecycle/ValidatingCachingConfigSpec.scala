package cromwell.engine.workflow.lifecycle

import com.typesafe.config.ConfigFactory
import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.{FlatSpec, Matchers}
import scala.util.{Failure, Success, Try}


class ValidatingCachingConfigSpec extends FlatSpec with Matchers with TableDrivenPropertyChecks {

    it should "run config tests" in {
      val cases = Table(
        ("config"                 , "exceptionMessage"                                                              ),
        ("enabled = not-a-boolean", "String: 1: enabled has type STRING rather than BOOLEAN"                        ),
        ("enabled = true"         , true                                                                            ),
        ("enabled = false"        , false                                                                           ),
        ("enabled = 1"            , "String: 1: enabled has type NUMBER rather than BOOLEAN"                        ),
        (""                       , "No configuration setting found for key 'enabled'"                              )
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


