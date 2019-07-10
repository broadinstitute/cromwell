package wom.util

import com.typesafe.config.ConfigFactory

object WomParseConfig {
  private val cfg = ConfigFactory.load().getConfig("wom-parse")
  val convertNestedScatterToSubworkflow = cfg.getBoolean("convert-nested-scatter-to-subworkflow")
}
