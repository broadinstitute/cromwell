package lenthall.util

import java.net.{MalformedURLException, URL}

import com.typesafe.config.{Config, ConfigValue}
import org.slf4j.LoggerFactory

import scala.collection.JavaConversions._
import scala.reflect._
import scalaz.Scalaz._
import scalaz._

/**
 * Enhances a Typesafe Config with the Scalaz validation utilities.
 *
 * If used in your project, scalaz must also be added to the build.sbt via:
 *   "org.scalaz" %% "scalaz-core" % "7.1.4"
 */
object ValidatedConfig {
  val validationLogger = LoggerFactory.getLogger("ConfigurationValidation")

  implicit class EnhancedValidatedConfig(val config: Config) extends AnyVal {
    def keys: Set[String] = {
      config.entrySet().toSet map { entry: java.util.Map.Entry[String, ConfigValue] => entry.getKey }
    }

    /**
     * For keys that are in the configuration but not in the reference keySet, log a warning.
     */
    def warnNotRecognized(keySet: Set[String], context: String) = {
      val unrecognizedKeys = keys.diff(keySet)
      if (unrecognizedKeys.nonEmpty) {
        validationLogger.warn(s"Unrecognized configuration key(s) for $context: ${unrecognizedKeys.mkString(", ")}")
      }
    }

    /**
     * Validates that the value for this key is a well formed URL.
     */
    def validateURL(key: String): ValidationNel[String, URL] = {
      if (config.hasPath(key)) {
        key.validateAny[URL, MalformedURLException] { url =>
          new URL(config.getString(url))
        }
      } else keyFailNel(key)
    }

    def validateString(key: String): ValidationNel[String, String] = {
      if (config.hasPath(key)) config.getString(key).successNel else keyFailNel(key)
    }

    def validateBoolean(key: String): ValidationNel[String, Boolean] = {
      if (config.hasPath(key)) config.getBoolean(key).successNel else keyFailNel(key)
    }

    def validateInt(key: String): ValidationNel[String, Int] = {
      if (config.hasPath(key)) config.getInt(key).successNel else keyFailNel(key)
    }

    def validateLong(key: String): ValidationNel[String, Long] = {
      if (config.hasPath(key)) config.getLong(key).successNel else keyFailNel(key)
    }

    def validateDouble(key: String): ValidationNel[String, Double] = {
      if (config.hasPath(key)) config.getDouble(key).successNel else keyFailNel(key)
    }

    @inline
    private def keyFailNel(key: String) = s"Could not find key: $key".failureNel
  }

  implicit class EnhancedValidation[I <: AnyRef](val value: I) extends AnyVal {
    /**
     * Validates this value by applying validationFunction to it and returning a Validation:
     * Returns successNel upon success.
     * If an exception is thrown and is a subtype of E, return failureNel with the exception message.
     * @param validationFunction function that should throw an exception if this value is found not to be valid
     * @tparam O return type of validationFunction
     * @tparam E Restricts the subtype of Exception that should be caught during validation
     */
    def validateAny[O, E <: Exception : ClassTag](validationFunction: I => O): ValidationNel[String, O] = try {
      validationFunction(value).successNel
    } catch {
      case e if classTag[E].runtimeClass.isInstance(e) => e.getMessage.failureNel
    }
  }

}
