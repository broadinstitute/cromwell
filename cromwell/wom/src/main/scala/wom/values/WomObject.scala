package wom.values

import wom.TsvSerializable
import wom.types._
import wom.util.FileUtil

import scala.util.{Failure, Success, Try}

trait WomObjectLike {
  def value: Map[String, WomValue]
  def copyWith(values: Map[String, WomValue]): WomValue with WomObjectLike
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
          attributes map { obj.value(_).valueString } mkString "\t"
        } mkString(start = "", sep = "\n", end = "\n")

        Success(s"$attributesLine\n$valuesLines")
      case _ => Failure(new UnsupportedOperationException("Could not serialize array: Objects in the array have different attributes."))
    }
  }

}

case class WomObject(value: Map[String, WomValue]) extends WomValue with WomObjectLike with TsvSerializable {
  val womType = WomObjectType

  override def toWomString: String =
    "object {" + value.map {case (k, v) => s"$k: ${v.toWomString}"}.mkString(", ") + "}"

  lazy val orderedAttributes = value.keySet.toSeq
  lazy val orderedValues = orderedAttributes map { value(_) }

  def tsvSerialize: Try[String] = Try {
    val keysLine = orderedAttributes.mkString(start = "", sep = "\t", end = "\n")
    val values = orderedValues map {
      case v if v.isInstanceOf[WomPrimitive] => v.valueString
      case _ => throw new UnsupportedOperationException("Can only TSV serialize an Object with Primitive values.")
    }

    keysLine + values.mkString(start = "", sep = "\t", end = "\n")
  }

  override def copyWith(values: Map[String, WomValue]) = copy(values)
}
