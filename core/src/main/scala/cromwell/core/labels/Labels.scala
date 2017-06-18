package cromwell.core.labels

import cats.data.Validated._
import cats.instances.vector._
import cats.syntax.traverse._
import lenthall.validation.ErrorOr.ErrorOr

import scala.collection.JavaConverters._

case class Labels(value: Vector[Label]) {

  def asMap = (value map { label => label.key -> label.value }).toMap

  def asJavaMap = asMap.asJava

  def ++(that: Labels) = Labels(value ++ that.value)
}

object Labels {
  def apply(values: (String, String)*): Labels = {
    val kvps: Seq[(String, String)] = values.toSeq
    Labels((kvps map { case (k, v) => Label.safeLabel(k, v) }).to[Vector])
  }

  def validateMapOfLabels(labels: Map[String, String]): ErrorOr[Labels] = {
    val vectorOfValidatedLabel: Vector[ErrorOr[Label]] = labels.toVector map { case (k, v) => Label.validateLabel(k, v) }

    vectorOfValidatedLabel.sequence[ErrorOr, Label] map { validatedVectorOfLabel => Labels(validatedVectorOfLabel) }
  }

  def empty = Labels(Vector.empty)
}
