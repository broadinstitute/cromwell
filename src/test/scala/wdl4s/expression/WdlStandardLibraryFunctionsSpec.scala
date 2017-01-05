package wdl4s.expression

import org.apache.commons.lang3.NotImplementedException
import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.{FlatSpec, Matchers}
import wdl4s.types.{WdlArrayType, WdlIntegerType, WdlOptionalType}
import wdl4s.values.{WdlArray, WdlFile, WdlFloat, WdlInteger, WdlOptionalValue, WdlString, WdlValue}

import scala.util.{Success, Try}


class WdlStandardLibraryFunctionsSpec extends FlatSpec with Matchers {

  import TableDrivenPropertyChecks._

  behavior of "WdlStandardLibraryFunctions' built-ins"

  val selectionTable = Table(
    ("input", "select_first", "select_all"),
    (List(None, Some(4), Some(5)), Some(4), List(4, 5)),
    (List(Some(3), None, Some(5)), Some(3), List(3, 5)),
    (List(None, None, None), None, List.empty[Int])
  )

  selectionTable foreach { case (input, select_first_ifAppropriate, select_all) =>

    val functionInput = Seq(Success(mkWdlArray(input)))
    if (select_first_ifAppropriate.isDefined) {
      val select_first = select_first_ifAppropriate.get
      it should s"select '$select_first' as the select_first value in [${input.mkString(", ")}]" in {
        val expectedSelectFirstOutput = Success(WdlInteger(select_first))
        TestableFunctions.select_first(functionInput) should be(expectedSelectFirstOutput)
      }
    }

    it should s"select [${select_all.mkString(", ")}] as the select_all value in [${input.mkString(", ")}]" in {
      val expectedSelectAllOutput = Success(WdlArray(WdlArrayType(WdlIntegerType), select_all.map(WdlInteger(_))))
      TestableFunctions.select_all(functionInput) should be(expectedSelectAllOutput)
    }
  }

  def mkWdlArray(ints: List[Option[Int]]): WdlArray = {
    val values = ints map {
      case Some(i) => WdlOptionalValue(WdlInteger(i))
      case None => WdlOptionalValue(WdlIntegerType, None)
    }
    WdlArray(WdlArrayType(WdlOptionalType(WdlIntegerType)), values)
  }
}

case object TestableFunctions extends WdlStandardLibraryFunctions {
  // No need to test the ones that are overridden anyway:
  // TODO: Can replace with "OnlyPureFunctions when that branch merges..."
  override def readFile(path: String): String = throw new NotImplementedException("")
  override def range(params: Seq[Try[WdlValue]]): Try[WdlArray] = throw new NotImplementedException("")
  override def read_json(params: Seq[Try[WdlValue]]): Try[WdlValue] = throw new NotImplementedException("")
  override def write_json(params: Seq[Try[WdlValue]]): Try[WdlFile] = throw new NotImplementedException("")
  override def sub(params: Seq[Try[WdlValue]]): Try[WdlString] = throw new NotImplementedException("")
  override def size(params: Seq[Try[WdlValue]]): Try[WdlFloat] = throw new NotImplementedException("")
  override def length(params: Seq[Try[WdlValue]]): Try[WdlInteger] = throw new NotImplementedException("")
  override def transpose(params: Seq[Try[WdlValue]]): Try[WdlArray] = throw new NotImplementedException("")
  override def write_tsv(params: Seq[Try[WdlValue]]): Try[WdlFile] = throw new NotImplementedException("")
  override def stdout(params: Seq[Try[WdlValue]]): Try[WdlFile] = throw new NotImplementedException("")
  override def glob(path: String, pattern: String): Seq[String] = throw new NotImplementedException("")
  override def stderr(params: Seq[Try[WdlValue]]): Try[WdlFile] = throw new NotImplementedException("")
  override def writeTempFile(path: String, prefix: String, suffix: String, content: String): String = throw new NotImplementedException("")
}

case object TestableFunctionTypes extends WdlStandardLibraryFunctionsType
