package lenthall.config

import com.typesafe.config.Config

object ScalaConfig {

  implicit class EnhancedScalaConfig(val config: Config) extends AnyVal {
    def getConfigOption(key: String): Option[Config] = {
      if (config.hasPath(key)) Option(config.getConfig(key)) else None
    }

    def getStringOption(key: String): Option[String] = {
      if (config.hasPath(key)) Option(config.getString(key)) else None
    }

    def getBooleanOption(key: String): Option[Boolean] = {
      if (config.hasPath(key)) Option(config.getBoolean(key)) else None
    }

    def getIntOption(key: String): Option[Int] = {
      if (config.hasPath(key)) Option(config.getInt(key)) else None
    }

    def getLongOption(key: String): Option[Long] = {
      if (config.hasPath(key)) Option(config.getLong(key)) else None
    }

    def getDoubleOption(key: String): Option[Double] = {
      if (config.hasPath(key)) Option(config.getDouble(key)) else None
    }

    def getStringOr(key: String, default: => String = ""): String = {
      getStringOption(key) getOrElse default
    }

    def getBooleanOr(key: String, default: => Boolean = false): Boolean = {
      getBooleanOption(key) getOrElse default
    }

    def getIntOr(key: String, default: => Int = 0): Int = {
      getIntOption(key) getOrElse default
    }

    def getLongOr(key: String, default: => Long = 0L): Long = {
      getLongOption(key) getOrElse default
    }

    def getDoubleOr(key: String, default: => Double = 0.0): Double = {
      getDoubleOption(key) getOrElse default
    }
  }

}
