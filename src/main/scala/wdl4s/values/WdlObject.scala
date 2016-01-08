package wdl4s.values

import wdl4s.types._
import wdl4s.{Call, TsvSerializable}
import wdl4s.util.FileUtil
import scala.util.{Failure, Success, Try}

trait WdlObjectLike {
  def value: Map[String, WdlValue]
}

object WdlObject {

  def coerceObject(m: Map[String, String]): WdlObject = {
    val coerced = WdlMap.coerceMap(m, WdlMapType(WdlStringType, WdlAnyType)).value map {
      case (k, v) => k.valueString -> v
    }

    WdlObject(coerced)
  }

  def fromTsv(tsv: String): Try[Array[WdlObject]] = {
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
  def tsvSerializeArray(input: Seq[WdlObject]): Try[String] = {

    /**
     * Validates that all objects have the same attributes.
     */
    def attributesMatch(objects: Seq[WdlObject]): Boolean = {
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
        } mkString "\n"

        Success(s"$attributesLine\n$valuesLines")
      case _ => Failure(new UnsupportedOperationException("Could not serialize array: Objects in the array have different attributes."))
    }
  }

}

case class WdlObject(value: Map[String, WdlValue]) extends WdlValue with WdlObjectLike with TsvSerializable {
  val wdlType = WdlObjectType

  override def toWdlString: String =
    "object {" + value.map {case (k, v) => s"$k: ${v.toWdlString}"}.mkString(", ") + "}"

  lazy val orderedAttributes = value.keySet.toSeq
  lazy val orderedValues = orderedAttributes map { value(_) }

  def tsvSerialize: Try[String] = Try {
    val keysLine = orderedAttributes.mkString("\t")
    val values = orderedValues map {
      case v if v.isInstanceOf[WdlPrimitive] => v.valueString
      case _ => throw new UnsupportedOperationException("Can only TSV serialize an Object with Primitive values.")
    }

    s"$keysLine\n${values.mkString("\t")}"
  }
}

case class WdlCallOutputsObject(call: Call, outputs: Map[String, WdlValue]) extends WdlValue with WdlObjectLike {
  val wdlType = WdlCallOutputsObjectType(call)
  val value = outputs
}
