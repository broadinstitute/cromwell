package cromwell.core.labels

import cats.data.Validated._
import cats.instances.vector._
import cats.syntax.traverse._
import lenthall.validation.ErrorOr
import lenthall.validation.ErrorOr.ErrorOr

import scala.collection.JavaConverters._

case class Labels(value: Vector[Label]) {

  def asTuple: Vector[(String, String)] = value.map(label => label.key -> label.value)

  def asMap: Map[String, String] = asTuple.toMap

  def asJavaMap = asMap.asJava

  def ++(that: Labels) = Labels(value ++ that.value)
}

object Labels {
  def apply(values: (String, String)*): Labels = {
    val kvps: Seq[(String, String)] = values.toSeq
    Labels((kvps map { case (k, v) => Label(k, v) }).to[Vector])
  }

  //This method should only ever be used by the JES backend as this manipulation
  //to safe google label names is required to guard against JES job submission failure.
  def googleLabels(values: (String, String)*): Labels = {

    def safeGoogleLabel(key: String, value: String): Label = {
      Label(Label.safeGoogleName(key), Label.safeGoogleName(value, emptyAllowed = true))
    }

    val kvps: Seq[(String, String)] = values.toSeq
    Labels((kvps map { case(k, v) => safeGoogleLabel(k, v) } ).to[Vector])
  }

  def validateMapOfLabels(labels: Map[String, String]): ErrorOr[Labels] = {
    labels.toVector traverse { Label.validateLabel _ }.tupled map Labels.apply
  }

  def empty = Labels(Vector.empty)
}
