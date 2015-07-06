package cromwell.engine.backend.jes

import cromwell.binding.values.{WdlInteger, WdlString, WdlGcsObject}
import org.scalatest.{Matchers, FlatSpec}

import scala.util.{Try, Success}

/**
 * Created by chrisl on 7/7/15.
 */
class JesEngineFunctionsSpec extends FlatSpec with Matchers{

  final val SECRETS_FILE_LOCATION = "/Users/chrisl/client_secrets.json"

  final val BUCKET_NAME = "chrisl-dsde-dev"

  final val INT_FILE = "intfile"
  final val INT_FILE_VALUE = 2012

  final val STRING_FILE = "BobLoblawsLawBlog"
  final val STRING_FILE_CONTENTS = """Day 1:
                                     |Law!
                                     |
                                     |Day 2:
                                     |Law!
                                     |
                                     |Day 3:
                                     |Law!
                                     |
                                     |Day 4:
                                     |Law!
                                     |
                                     |Day 5:
                                     |Law!
                                     |
                                     |Day 6:
                                     |Law!
                                     |
                                     |Day 7:
                                     |Rest.
                                     |""".stripMargin

  "JES Engine Functions" should " read strings correctly" in {
    val readString = new JesEngineFunctions(SECRETS_FILE_LOCATION).getFunction("read_string")
    val gcsPathTry: Try[WdlGcsObject] = Success(WdlGcsObject(BUCKET_NAME, STRING_FILE))
    readString(Seq(gcsPathTry)) shouldEqual Success(WdlString(STRING_FILE_CONTENTS))
  }

  "JES Engine Functions" should " read ints correctly" in {
    val readString = new JesEngineFunctions(SECRETS_FILE_LOCATION).getFunction("read_int")
    val gcsPathTry: Try[WdlGcsObject] = Success(WdlGcsObject(BUCKET_NAME, INT_FILE))
    readString(Seq(gcsPathTry)) shouldEqual Success(WdlInteger(INT_FILE_VALUE))
  }
}
