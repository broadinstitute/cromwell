package cromwell.services.loadcontroller.impl

import cromwell.services.loadcontroller.LoadMetric

case class MetadataQueueMetric(threshold: Int) extends LoadMetric {
  val name = "metadataQueueSize"
  val highLoadThreshold: Int = threshold
}
