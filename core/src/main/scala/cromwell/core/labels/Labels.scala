package cromwell.core.labels

import scala.collection.JavaConverters._

case class Labels(value: Vector[Label]) {

  def asJesLabels = (value map { label => label.key -> label.value }).toMap.asJava

  def ++(that: Labels) = Labels(value ++ that.value)
}

object Labels {
  def apply(values: (String, String)*): Labels = {
    val kvps: Seq[(String, String)] = values.toSeq
    Labels((kvps map { case (k, v) => Label.safeLabel(k, v) }).to[Vector])
  }

  def empty = Labels(Vector.empty)
}
