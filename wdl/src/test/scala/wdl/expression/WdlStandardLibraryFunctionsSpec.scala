package wdl.expression

import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.{FlatSpec, Matchers}
import wom.types.{WomArrayType, WomIntegerType, WomOptionalType}
import wom.values._

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
        val expectedSelectFirstOutput = Success(WomInteger(select_first))
        TestableFunctions.select_first(functionInput) should be(expectedSelectFirstOutput)
      }
    }

    it should s"select [${select_all.mkString(", ")}] as the select_all value in [${input.mkString(", ")}]" in {
      val expectedSelectAllOutput = Success(WomArray(WomArrayType(WomIntegerType), select_all.map(WomInteger)))
      TestableFunctions.select_all(functionInput) should be(expectedSelectAllOutput)
    }
  }

  def mkWdlArray(ints: List[Option[Int]]): WomArray = {
    val values = ints map {
      case Some(i) => WomOptionalValue(WomInteger(i))
      case None => WomOptionalValue(WomIntegerType, None)
    }
    WomArray(WomArrayType(WomOptionalType(WomIntegerType)), values)
  }
}

case object TestableFunctions extends WdlStandardLibraryFunctions {
  // No need to test the ones that are overridden anyway:
  // TODO: Can replace with "OnlyPureFunctions when that branch merges..."
  override def readFile(path: String): String = throw new NotImplementedError
  override def writeFile(path: String, content: String): Try[WomFile] = throw new NotImplementedError
  override def range(params: Seq[Try[WomValue]]): Try[WomArray] = throw new NotImplementedError
  override def read_json(params: Seq[Try[WomValue]]): Try[WomValue] = throw new NotImplementedError
  override def write_json(params: Seq[Try[WomValue]]): Try[WomFile] = throw new NotImplementedError
  override def sub(params: Seq[Try[WomValue]]): Try[WomString] = throw new NotImplementedError
  override def size(params: Seq[Try[WomValue]]): Try[WomFloat] = throw new NotImplementedError
  override def length(params: Seq[Try[WomValue]]): Try[WomInteger] = throw new NotImplementedError
  override def transpose(params: Seq[Try[WomValue]]): Try[WomArray] = throw new NotImplementedError
  override def write_tsv(params: Seq[Try[WomValue]]): Try[WomFile] = throw new NotImplementedError
  override def stdout(params: Seq[Try[WomValue]]): Try[WomFile] = throw new NotImplementedError
  override def globHelper(pattern: String): Seq[String] = throw new NotImplementedError
  override def stderr(params: Seq[Try[WomValue]]): Try[WomFile] = throw new NotImplementedError
}

case object TestableFunctionTypes extends WdlStandardLibraryFunctionsType
