package cromwell.util

import java.net.URL

import com.typesafe.config.{ConfigValue, Config, ConfigException}
import org.slf4j.{LoggerFactory, Logger}

import scala.collection.JavaConversions._
import scala.util.Try

object ConfigUtil {

  val validationLogger = LoggerFactory.getLogger("ConfigurationValidation")
  
  class ConfigValidationException(context: String, validationException: ConfigException.ValidationFailed)
    extends ConfigException.ValidationFailed(validationException.problems()) {
    override def getMessage: String = {
      val problems = validationException.problems().map(_.problem()).mkString(", ")
      s"$context: Configuration Validation failed: $problems"
    }
  }

  implicit class EnhancedConfig(val config: Config) extends AnyVal {
    
    def getStringOption(key: String): Option[String] = {
      Try(config.getString(key)) match {
        case scala.util.Success(value) => Option(value)
        case scala.util.Failure(e: ConfigException.Missing) => None
        case scala.util.Failure(e) => throw e
      }
    }
    
    def getURL(key: String): URL = {
      new URL(config.getString(key))
    }

    /**
     * Checks that the configuration is valid compared to refConfig (i.e. at least all same keys with same type).
     * For keys that are in the configuration but not in the reference configuration, log a warning.
     * @param optionalKeys keys that should not be considered unrecognized if they are present but are not in the reference configuration
     */
    def checkValidAndWarnNotRecognized(refConfig: Config, context: String, optionalKeys: String*) = {
      try {
        config.checkValid(refConfig)
      } catch {
        case e: ConfigException.ValidationFailed => throw new ConfigValidationException(context, e)
        case t: Throwable => throw t
      }
      warnNotRecognized(config, refConfig, context, optionalKeys)
    }

  }
  
  def warnNotRecognized(config: Config, refConf: Config, context: String, optionalKeys: Seq[String]) = {
    val refKeys = refConf.entrySet().toSet map { v: java.util.Map.Entry[String, ConfigValue] => v.getKey }
    val confKeys = config.entrySet().toSet map { v: java.util.Map.Entry[String, ConfigValue] => v.getKey }

    confKeys.diff(refKeys ++ optionalKeys) match {
      case warnings if warnings.nonEmpty => validationLogger.warn(s"Unrecognized configuration key for $context: ${warnings.mkString(", ")}")
      case _ =>
    }
  }
}
