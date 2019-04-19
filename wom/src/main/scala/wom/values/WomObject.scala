package wom.values

import cats.Applicative
import cats.instances.list._
import cats.syntax.traverse._
import cats.syntax.functor._
import common.Checked
import common.validation.ErrorOr.ErrorOr
import common.validation.IOChecked.IOChecked
import wom.TsvSerializable
import wom.expression.IoFunctionSet
import wom.types._
import wom.util.FileUtil

import scala.language.higherKinds
import scala.util.{Failure, Success, Try}

trait WomObjectLike extends WomValue {
  def values: Map[String, WomValue]
  def copyWith(values: Map[String, WomValue]): WomValue with WomObjectLike
  override def toWomString: String = {
    val typed = womObjectTypeLike match {
      case _: WomCompositeType => "typed "
      case _ => ""
    }
    s"${typed}object {" + values.map { case (k, v) => s"$k: ${v.toWomString}" }.mkString(", ") + "}"
  }
  def womObjectTypeLike: WomObjectTypeLike

  override def collectAsSeq[T <: WomValue](filterFn: PartialFunction[WomValue, T]): Seq[T] = {
    values.values.toSeq flatMap { _.collectAsSeq(filterFn) }
  }

  def traverse[R <: WomValue, G[_]](f: WomValue => G[R])(implicit applicative: Applicative[G]): G[WomObjectLike] = {
    if (values.isEmpty) applicative.pure(this)
    else {
      val traverseFunction: (String, WomValue) => G[(String, R)] = {
        case (key, value) => applicative.map(f(value)) { key -> _ }
      }
      applicative.map(values.toList.traverse[G, (String, R)](traverseFunction.tupled)) { mapped =>
        copyWith(mapped.toMap)
      }
    }
  }

  override def initialize(ioFunctionSet: IoFunctionSet): IOChecked[WomValue] = traverse(_.initialize(ioFunctionSet)).widen
}

object WomObject {

  def coerceObject(m: Map[String, String]): WomObject = {
    val coerced = WomMap.coerceMap(m, WomMapType(WomStringType, WomAnyType)).value map {
      case (k, v) => k.valueString -> v
    }

    WomObject(coerced)
  }

  def fromTsv(tsv: String): Try[Array[WomObject]] = {
    FileUtil.parseTsv(tsv) match {
      case Success(table) if table.isEmpty => Failure(new UnsupportedOperationException("TSV file was empty or could not be parsed."))
      case Success(table) if table.length < 2 => Failure(new UnsupportedOperationException("TSV must be 2 rows (or more) to convert to an Object (Array[Object])"))
      case Success(table) => Try {
        table.tail map { line => coerceObject((table.head zip line).toMap) }
      }
      case Failure(e) => Failure(e)
    }
  }

  //TODO: Try to stream this out to avoid memory overhead
  def tsvSerializeArray(input: Seq[WomObject]): Try[String] = {

    // Validates that all objects have the same attributes.
    def attributesMatch(objects: Seq[WomObject]): Boolean = {
      val attributesSet = objects map { _.orderedAttributes }
      val intersection = attributesSet reduce (_.intersect(_))
      attributesSet forall { _ == intersection }
    }

    input match {
      case Nil => Failure(new UnsupportedOperationException("Cannot write empty objects array."))
      case objects if attributesMatch(objects) =>
        /* Note: this is arbitrary as it takes the first object in the array as a reference for the attributes.
         * It has no impact though as we first made sure that all objects have the same attributes.
         */
        val attributes = objects.head.orderedAttributes
        val attributesLine = attributes.mkString("\t")
        val valuesLines = objects map { obj =>
          attributes map { obj.values(_).valueString } mkString "\t"
        } mkString(start = "", sep = "\n", end = "\n")

        Success(s"$attributesLine\n$valuesLines")
      case _ => Failure(new UnsupportedOperationException("Could not serialize array: Objects in the array have different attributes."))
    }
  }
  
  def apply(values: Map[String, WomValue]) = new WomObject(values, WomObjectType)
  
  def withTypeUnsafe(values: Map[String, Any], objectTypeLike: WomObjectTypeLike): WomObject = {
    import common.validation.Validation._
    withTypeErrorOr(values, objectTypeLike).toTry.get
  }

  def withTypeErrorOr(values: Map[String, Any], objectTypeLike: WomObjectTypeLike): ErrorOr[WomObject] = {
    objectTypeLike.validateAndCoerceValues(values).map(new WomObject(_, objectTypeLike))
  }

  def withTypeChecked(values: Map[String, Any], objectTypeLike: WomObjectTypeLike): Checked[WomObject] = {
    withTypeErrorOr(values, objectTypeLike).toEither
  }

}

case class WomObject private[WomObject] (values: Map[String, WomValue], womType: WomObjectTypeLike) extends WomObjectLike with TsvSerializable {
  lazy val orderedAttributes = values.keySet.toSeq
  lazy val orderedValues = orderedAttributes map { values(_) }
  lazy val womObjectTypeLike = womType
  
  def tsvSerialize: Try[String] = Try {
    val keysLine = orderedAttributes.mkString(start = "", sep = "\t", end = "\n")
    val values = orderedValues map {
      case v: WomPrimitive => v.valueString
      case _ => throw new UnsupportedOperationException("Can only TSV serialize an Object with Primitive values.")
    }

    keysLine + values.mkString(start = "", sep = "\t", end = "\n")
  }

  override def copyWith(values: Map[String, WomValue]) = copy(values)
}
