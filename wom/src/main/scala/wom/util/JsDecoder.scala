package wom.util

import cats.instances.list._
import cats.syntax.traverse._
import cats.syntax.validated._
import common.validation.ErrorOr._
import org.mozilla.javascript.{NativeArray, NativeObject}
import wom.types._
import wom.util.JsUtil.{JsArray, JsObject}
import wom.values.{WomArray, WomBoolean, WomFloat, WomInteger, WomObject, WomOptionalValue, WomString, WomValue}

import scala.collection.JavaConverters._

/**
  * Converts a javascript value into a WomValue.
  *
  * Uses a lot of runtime pattern matching / partial functions.
  *
  * Perhaps a future version could use compile time checks with shapeless / polys.
  */
class JsDecoder {

  /**
    * Decodes a generic javascript type.
    */
  def decode(value: AnyRef): ErrorOr[WomValue] = {
    value match {

      case arr: NativeArray => arr.toArray.map(_.asInstanceOf[NativeObject]).toList.traverse[ErrorOr, WomValue](decode).map(WomArray)
      case obj: NativeObject if (obj.getClassName.contains("File")) =>
//        obj.entrySet().asScala.map {
//        case (key, value) =>
//      }
      case other =>
        val y = other
        println(y)
        ???
    }
  }
}

