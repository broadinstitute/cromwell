package cromwell.core.io

import scala.concurrent.duration.FiniteDuration

case class Throttle(elements: Int, per: FiniteDuration, maximumBurst: Int)
