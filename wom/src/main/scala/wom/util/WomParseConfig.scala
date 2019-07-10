package wom.util

import com.typesafe.config.ConfigFactory

object WomParseConfig {
  private val cfg = ConfigFactory.load().getConfig("womParse")
  val convertNestedScatterToSubworkflow = cfg.getBoolean("convert-nested-scatter-to-subworkflow")
}
