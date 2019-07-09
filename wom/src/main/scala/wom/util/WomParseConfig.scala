package wom.util

import com.typesafe.config.ConfigFactory

object WomParseConfig {
  private val cfg = ConfigFactory.load().getConfig("womParse")
  val innerOuterScatter = cfg.getBoolean("inner-outer-scatter")
}
