package cromwell.engine.backend.jes

// FIXME: See DSDEEPB-828, leaving this here as an example for when that ticket gets picked up

///**
// * Specification for the JES Engine Functions
// */
//class JesEngineFunctionsSpec extends FlatSpec with Matchers{
//
//  val secretsFileLocation: Path = Paths.get("/Users/chrisl/client_secrets.json")
//
//  final val BUCKET_NAME = "chrisl-dsde-dev"
//
//  final val INT_FILE = "intfile"
//  final val INT_FILE_VALUE = 2012
//
//  final val STRING_FILE = "BobLoblawsLawBlog"
//  final val STRING_FILE_CONTENTS = """Day 1:
//                                     |Law!
//                                     |
//                                     |Day 2:
//                                     |Law!
//                                     |
//                                     |Day 3:
//                                     |Law!
//                                     |
//                                     |Day 4:
//                                     |Law!
//                                     |
//                                     |Day 5:
//                                     |Law!
//                                     |
//                                     |Day 6:
//                                     |Law!
//                                     |
//                                     |Day 7:
//                                     |Rest.
//                                     |""".stripMargin
//
//  "JES Engine Functions" should "read strings correctly" in {
//    val readString = JesEngineFunctions(secretsFileLocation, "gs://a/a").getFunction("read_string")
//    val gcsPathTry: Try[WdlFile] = Success(WdlFile(GoogleCloudStoragePath(BUCKET_NAME, STRING_FILE).toString))
//    readString(Seq(gcsPathTry)) shouldEqual Success(WdlString(STRING_FILE_CONTENTS))
//  }
//
//  "JES Engine Functions" should " read ints correctly" in {
//    val readString = JesEngineFunctions(secretsFileLocation, "gs://a/a").getFunction("read_int")
//    val gcsPathTry: Try[WdlFile] = Success(WdlFile(GoogleCloudStoragePath(BUCKET_NAME, INT_FILE).toString))
//    readString(Seq(gcsPathTry)) shouldEqual Success(WdlInteger(INT_FILE_VALUE))
//  }
//
//  "JES Engine Functions" should "be able to generate valid GCS paths from random strings" in {
//    val randomString = "whoopwhoop!!"
//    val gcsPath = JesEngineFunctions(secretsFileLocation, "gs://a/a").gcsPathFromAnyString(randomString)
//    assert(gcsPath.bucket equals "a")
//    assert(gcsPath.objectName equals "a/MkB+wg3Z1ZS/m9trfL7qNw==")
//  }
//}
