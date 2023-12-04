package cromwell.core.labels

import cats.data.Validated._
import cats.instances.vector._
import cats.syntax.traverse._
import common.validation.ErrorOr.ErrorOr

import scala.jdk.CollectionConverters._

case class Labels(value: Vector[Label]) {

  def asTuple: Vector[(String, String)] = value.map(label => label.key -> label.value)

  def asMap: Map[String, String] = asTuple.toMap

  def asJavaMap = asMap.asJava

  def ++(that: Labels) = Labels(value ++ that.value)
}

object Labels {
  def apply(values: (String, String)*): Labels =
    Labels(values.toVector map (Label.apply _).tupled)

  def validateMapOfLabels(labels: Map[String, String]): ErrorOr[Labels] =
    labels.toVector traverse { Label.validateLabel _ }.tupled map Labels.apply

  def empty = Labels(Vector.empty)
}
