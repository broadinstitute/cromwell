package cromwell.api

import akka.actor.ActorSystem
import akka.http.scaladsl.model.{ContentType, HttpEntity}
import akka.stream.ActorMaterializer
import akka.stream.scaladsl.Sink
import better.files.File
import cromwell.api.model.{Label, WorkflowBatchSubmission, WorkflowSingleSubmission}
import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.{AsyncFlatSpec, BeforeAndAfterAll, Matchers}
import spray.json.JsonParser.ParsingException

class CromwellClientSpec extends AsyncFlatSpec with BeforeAndAfterAll with Matchers with TableDrivenPropertyChecks {
  behavior of "CromwellClient"

  implicit val system = ActorSystem("CromwellClientSpec")
  implicit val materializer = ActorMaterializer()

  private val tempFile: File = File.newTemporaryFile("cromwell_client_spec.", ".tmp").write("hello")

  override protected def afterAll(): Unit = {
    system.terminate()
    tempFile.delete(swallowIOExceptions = true)
    super.afterAll()
  }

  val okRefreshTokenTests = Table(
    ("description", "optionsOption", "refreshTokenOption", "expected"),
    ("ignore bad json when refresh token not provided", Option("{"), None, Option("{")),
    ("ignore bad json when refresh token provided but not used", Option("{"), Option("myToken"), Option("{")),
    ("not format json when refresh token key not found", Option("{   }"), Option("myToken"), Option("{   }")),
    ("replace token when found", Option("""{"refresh_token" : "replace_me"}"""), Option("myToken"),
      Option("""{"refresh_token":"myToken"}""")),
  )

  forAll(okRefreshTokenTests) { (description, optionsOption, refreshTokenOption, expected) =>
    it should description in {
      val actual = CromwellClient.replaceJson(optionsOption, "refresh_token", refreshTokenOption)
      actual should be(expected)
      succeed
    }
  }

  it should "throw an exception when inserting a refresh token into bad json using the token" in {
    val optionsOption = Option("""{"refresh_token" : "replace_me"""")
    val refreshTokenOption = Option("myToken")
    val actual = intercept[ParsingException] {
      CromwellClient.replaceJson(optionsOption, "refresh_token", refreshTokenOption)
    }
    actual.summary should be("""Unexpected end-of-input at input index 31 (line 1, position 32), expected '}'""")
    actual.detail should be(
      """|
         |{"refresh_token" : "replace_me"
         |                               ^
         |""".stripMargin)
    succeed
  }

  val okRequestEntityTests = Table(
    ("description", "workflowSubmission", "expectedJsons", "expectedFiles"),

    ("submit a wdl",
      WorkflowSingleSubmission("wdl", None, None, None, None, None, None),
      Map("workflowSource" -> "wdl"),
      Map()
    ),

    ("batch submit a wdl",
      WorkflowBatchSubmission("wdl", None, None, List(), None, None, None),
      Map("workflowSource" -> "wdl", "workflowInputs" -> "[]"),
      Map()
    ),

    ("submit a wdl with data",
      WorkflowSingleSubmission(
        "wdl",
        Option("wfType"),
        Option("wfTypeVersion"),
        Option("inputsJson"),
        Option("optionsJson"),
        Option(List(Label("labelKey", "labelValue"))),
        Option(tempFile)
      ),
      Map(
        "workflowSource" -> "wdl",
        "workflowType" -> "wfType",
        "workflowTypeVersion" -> "wfTypeVersion",
        "workflowInputs" -> "inputsJson",
        "workflowOptions" -> "optionsJson",
        "customLabels" -> """{"labelKey":"labelValue"}"""
      ),
      Map("workflowDependencies" -> tempFile)
    ),

    ("batch submit a wdl with data",
      WorkflowBatchSubmission(
        "wdl",
        Option("wfType"),
        Option("wfTypeVersion"),
        List("inputsJson1", "inputsJson2"),
        Option("optionsJson"),
        Option(List(Label("labelKey", "labelValue"))),
        Option(tempFile)
      ),
      Map(
        "workflowSource" -> "wdl",
        "workflowType" -> "wfType",
        "workflowTypeVersion" -> "wfTypeVersion",
        "workflowInputs" -> "[inputsJson1,inputsJson2]",
        "workflowOptions" -> "optionsJson",
        "customLabels" -> """{"labelKey":"labelValue"}"""
      ),
      Map("workflowDependencies" -> tempFile)
    )
  )

  forAll(okRequestEntityTests) { (description, workflowSubmission, expectedJsons, expectedFiles) =>
    it should description in {
      val actual = CromwellClient.requestEntityForSubmit(workflowSubmission)
      actual should be(a[HttpEntity.Chunked])
      actual match {
        case HttpEntity.Chunked(contentType: ContentType.Binary, chunks) =>
          contentType.mediaType.isMultipart should be(true)
          val boundary = contentType.mediaType.params("boundary")

          val expectedJsonChunks = expectedJsons map {
            case (chunkKey, chunkValue) =>
              s"""|--$boundary
                  |Content-Type: application/json
                  |Content-Disposition: form-data; name="$chunkKey"
                  |
                  |$chunkValue
                  |""".stripMargin.replace("\n", "\r\n").trim
          }
          val expectedFileChunks = expectedFiles map {
            case (chunkKey, chunkFile) =>
              s"""|--$boundary
                  |Content-Type: application/zip
                  |Content-Disposition: form-data; filename="${chunkFile.name}"; name="$chunkKey"
                  |""".stripMargin.replace("\n", "\r\n").trim
          }
          val expectedFileContents = expectedFiles map {
            case (_, chunkFile) => chunkFile.contentAsString
          }
          val boundaryEnd = s"--$boundary--"

          val expectedChunks = Seq(boundaryEnd) ++ expectedJsonChunks ++ expectedFileChunks ++ expectedFileContents

          val futureChunks = chunks.map(_.data().utf8String).runWith(Sink.seq)
          futureChunks map { chunks =>
            val actualChunks = chunks.map(_.trim)
            actualChunks should contain theSameElementsAs expectedChunks
          } map { _ => succeed }
        case otherEntity =>
          fail(s"Unexpected entity $otherEntity")
      }
    }
  }
}
