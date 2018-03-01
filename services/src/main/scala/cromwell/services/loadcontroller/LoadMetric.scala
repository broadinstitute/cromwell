package cromwell.services.loadcontroller

import cromwell.services.loadcontroller.LoadControllerService._

trait LoadMetric {
  def name: String
  def highLoadThreshold: Int
  def loadLevel(value: Int) = if (value > highLoadThreshold) HighLoad else NormalLoad
}
