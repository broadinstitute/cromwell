package cromwell.api

import akka.actor.ActorSystem
import akka.http.scaladsl.model.{ContentType, HttpEntity}
import akka.stream.ActorMaterializer
import akka.stream.scaladsl.Sink
import better.files.File
import cromwell.api.model.{Label, WorkflowBatchSubmission, WorkflowId, WorkflowSingleSubmission}
import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.BeforeAndAfterAll
import org.scalatest.flatspec.AsyncFlatSpec
import org.scalatest.matchers.should.Matchers


class CromwellClientSpec extends AsyncFlatSpec with BeforeAndAfterAll with Matchers with TableDrivenPropertyChecks {
  behavior of "CromwellClient"

  implicit val system: ActorSystem = ActorSystem("CromwellClientSpec")
  implicit val materializer: ActorMaterializer = ActorMaterializer()

  private val tempFile: File = File.newTemporaryFile("cromwell_client_spec.", ".tmp").write("hello")

  override protected def afterAll(): Unit = {
    system.terminate()
    tempFile.delete(swallowIOExceptions = true)
    super.afterAll()
  }
  
  it should "build an uri with query arguments" in {
    val id = WorkflowId.randomId()
    val args = Option(
      Map(
        "key1" -> List("%v11%", "v12"),
        "key2" -> List("v2")
      )
    )
    CromwellClient.workflowSpecificGetEndpoint("http://submit", id, "endpoint", args).toString() shouldBe
      s"http://submit/$id/endpoint?key1=%25v11%25&key1=v12&key2=v2"
  }

  private val okRequestEntityTests = Table(
    ("description", "workflowSubmission", "expectedJsons", "expectedFiles"),

    ("submit a wdl",
      WorkflowSingleSubmission(Option("wdl"), None, None, None, None, None, None, None, None),
      Map("workflowSource" -> "wdl"),
      Map()
    ),

    ("batch submit a wdl",
      WorkflowBatchSubmission(Option("wdl"), None, None, None, None, List(), None, None, None),
      Map("workflowSource" -> "wdl", "workflowInputs" -> "[]"),
      Map()
    ),

    ("submit a wdl with data",
      WorkflowSingleSubmission(
        Option("wdl"),
        None,
        None,
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
        "labels" -> """{"labelKey":"labelValue"}"""
      ),
      Map("workflowDependencies" -> tempFile)
    ),

    ("submit a wdl using workflow url",
      WorkflowSingleSubmission(
        None,
        Option("https://link-to-url"),
        None,
        Option("wfType"),
        Option("wfTypeVersion"),
        Option("inputsJson"),
        Option("optionsJson"),
        Option(List(Label("labelKey", "labelValue"))),
        Option(tempFile)
      ),
      Map(
        "workflowUrl" -> "https://link-to-url",
        "workflowType" -> "wfType",
        "workflowTypeVersion" -> "wfTypeVersion",
        "workflowInputs" -> "inputsJson",
        "workflowOptions" -> "optionsJson",
        "labels" -> """{"labelKey":"labelValue"}"""
      ),
      Map("workflowDependencies" -> tempFile)
    ),

    ("batch submit a wdl with data",
      WorkflowBatchSubmission(
        Option("wdl"),
        None,
        None,
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
        "labels" -> """{"labelKey":"labelValue"}"""
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
