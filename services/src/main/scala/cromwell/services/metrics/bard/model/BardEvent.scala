package cromwell.services.metrics.bard.model

import scala.jdk.CollectionConverters._

trait BardEvent extends Product {

  private val baseProperties: Map[String, Any] = Map("event" -> eventName, "pushToMixpanel" -> false)

  def eventName: String

  def getProperties: java.util.Map[String, Any] =
    (this.productElementNames.zip(this.productIterator).toMap ++ baseProperties).asJava

}
