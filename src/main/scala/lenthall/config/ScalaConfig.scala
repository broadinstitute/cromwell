package lenthall.config

import com.typesafe.config.{Config, ConfigException, ConfigFactory}

import scala.util.{Failure, Success, Try}

object ScalaConfig {
  implicit class EnhancedScalaConfig(val config: Config) extends AnyVal {
    def getConfigOption(key: String): Option[Config] = getOption(key, config.getConfig)
    def getStringOption(key: String): Option[String] = getOption(key, config.getString)
    def getBooleanOption(key: String): Option[Boolean] = getOption(key, config.getBoolean)
    def getIntOption(key: String): Option[Int] = getOption(key, config.getInt)
    def getLongOption(key: String): Option[Long] = getOption(key, config.getLong)
    def getDoubleOption(key: String): Option[Double] = getOption(key, config.getDouble)
    def getBytesOption(key: String): Option[Long] = getOption(key, config.getBytes)
    def getConfigOr(key: String, default: => Config = empty): Config = getConfigOption(key) getOrElse default
    def getStringOr(key: String, default: => String = ""): String = getStringOption(key) getOrElse default
    def getBooleanOr(key: String, default: => Boolean = false): Boolean = getBooleanOption(key) getOrElse default
    def getIntOr(key: String, default: => Int = 0): Int = getIntOption(key) getOrElse default
    def getLongOr(key: String, default: => Long = 0L): Long = getLongOption(key) getOrElse default
    def getDoubleOr(key: String, default: => Double = 0.0): Double = getDoubleOption(key) getOrElse default
    def getBytesOr(key: String, default: => Long = 0L): Long = getBytesOption(key) getOrElse default
  }

  private[config] val empty = ConfigFactory.empty("empty")

  private[config] def getOption[T](key: String, f: String => T): Option[T] = {
    Try(f(key)) match {
      case Success(value) => Option(value)
      case Failure(e: ConfigException.Missing) => None
      case Failure(e) => throw e
    }
  }
}
