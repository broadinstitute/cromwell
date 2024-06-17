package cromwell.services.metrics.bard.model

trait BardEvent extends Product {

  private val baseProperties: Map[String, Any] = Map("event" -> eventName, "pushToMixpanel" -> false)

  def eventName: String

  def assembleScalaProperties: Map[String, Any] =
    this.productElementNames.zip(this.productIterator).toMap ++ baseProperties

  def getProperties: java.util.Map[String, Any]

}
