package com.readytalk.metrics

import org.slf4j.LoggerFactory
object CromwellStatsD {
  val logger = LoggerFactory.getLogger("StatsDLogger")
}

/**
  * Filters out unwanted timing metrics.
  * The package is funky on purpose because the StatsD constructor is package private so we need to be in this package to invoke it.
  * Taken from rawls.
  */
case class CromwellStatsD(hostname: String, port: Int) extends StatsD(hostname, port) {
  val MetricSuffixesToFilter = Set("max", "min", "p50", "p75", "p98", "p99", "p999", "mean_rate", "m5_rate", "m15_rate")
  
  override def send(name: String, value: String): Unit = {
    if (MetricSuffixesToFilter.exists(suffix => name.endsWith(suffix)))
      CromwellStatsD.logger.debug(s"Filtering metric with name [$name] and value [$value]")
    else super.send(name, value)
  }
}
