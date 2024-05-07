package org.broadinstitute.dsde.workbench.cromwell.consumer

import au.com.dius.pact.consumer.dsl.LambdaDsl.newJsonBody
import au.com.dius.pact.consumer.dsl._
import au.com.dius.pact.consumer.{ConsumerPactBuilder, PactTestExecutionContext}
import au.com.dius.pact.core.model.RequestResponsePact
import cats.effect.IO
import org.broadinstitute.dsde.workbench.cromwell.consumer.PactHelper._
import org.http4s.Uri
import org.http4s.blaze.client.BlazeClientBuilder
import org.http4s.client.Client
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import pact4s.scalatest.RequestResponsePactForger

import java.util.concurrent.Executors
import scala.concurrent.ExecutionContext

class DrsHubClientSpec extends AnyFlatSpec with Matchers with RequestResponsePactForger {
  val ec = ExecutionContext.fromExecutorService(Executors.newCachedThreadPool())
  implicit val cs = IO.contextShift(ec)
  /*
    Define the folder that the pact contracts get written to upon completion of this test suite.
   */
  override val pactTestExecutionContext: PactTestExecutionContext =
    new PactTestExecutionContext(
      "./target/pacts"
    )

  private val requestFields = List(
    "bucket",
    "accessUrl",
    "googleServiceAccount",
    "fileName",
    "hashes",
    "localizationPath",
    "bondProvider",
    "name",
    "size",
    "timeCreated",
    "timeUpdated",
    "gsUri",
    "contentType"
  )

  val filesize = 123L
  val timeCreated = "2021-03-04T20:00:00.000Z"
  val bucket = "fc-secure-1234567890"
  val filename = "my-file.bam"
  val bondProvider = "anvil"
  val fileHash = "a2317edbd2eb6cf6b0ee49cb81e3a556"
  val accessUrl = f"gs://${bucket}/${filename}"

  val drsResourceResponsePlaceholder: ResourceMetadata = ResourceMetadata(
    "application/octet-stream",
    filesize,
    timeCreated,
    timeCreated,
    None,
    None,
    None,
    None,
    None,
    Option(AccessUrl(accessUrl, List("Header", "Example"))),
    Map("md5" -> fileHash),
    None,
    Option(bondProvider)
  )

  val resourceMetadataResponseDsl: DslPart = newJsonBody { o =>
    o.stringType("contentType", "application/octet-stream")
    o.numberType("size", filesize)
    o.stringType("timeCreated", timeCreated)
    o.stringType("timeUpdated", timeCreated)
    o.nullValue("gsUri")
    o.nullValue("googleServiceAccount")
    o.nullValue("fileName")
    o.`object`("accessUrl",
               { a =>
                 a.stringType("url", accessUrl)
                 a.`array`("headers",
                           { h =>
                             h.stringType("Header")
                             h.stringType("Example")
                             ()
                           }
                 )
                 ()
               }
    )
    o.`object`("hashes",
               { o =>
                 o.stringType("md5", fileHash)
                 ()
               }
    )
    o.nullValue("localizationPath")
    o.stringType("bondProvider", bondProvider)
    ()
  }.build

  val fileId = "1234567890"

  val resourceRequestDsl = newJsonBody { o =>
    o.stringType("url", f"drs://test.theanvil.io/${fileId}")
    o.array("fields",
            { a =>
              requestFields.map(a.stringType)
              ()
            }
    )
    ()
  }.build

  val consumerPactBuilder: ConsumerPactBuilder = ConsumerPactBuilder
    .consumer("cromwell")

  val pactProvider: PactDslWithProvider = consumerPactBuilder
    .hasPactWith("drshub")

  var pactDslResponse: PactDslResponse = buildInteraction(
    pactProvider,
    state = "resolve Drs url",
    stateParams = Map[String, String](
      "fileId" -> fileId,
      "bucket" -> bucket,
      "filename" -> filename,
      "bondProvider" -> bondProvider,
      "fileHash" -> fileHash,
      "accessUrl" -> accessUrl,
      "fileSize" -> filesize.toString,
      "timeCreated" -> timeCreated
    ),
    uponReceiving = "Request to resolve drs url",
    method = "POST",
    path = "/api/v4/drs/resolve",
    requestHeaders = Seq("Content-Type" -> "application/json"),
    requestBody = resourceRequestDsl,
    status = 200,
    responseHeaders = Seq(),
    responseBody = resourceMetadataResponseDsl
  )

  pactDslResponse = buildInteraction(
    pactDslResponse,
    state = "Drshub is ok",
    uponReceiving = "Request for drshub api status",
    method = "GET",
    path = "/status",
    requestHeaders = Seq(),
    status = 200,
    responseHeaders = Seq()
  )

  override val pact: RequestResponsePact = pactDslResponse.toPact

  val client: Client[IO] =
    BlazeClientBuilder[IO](ExecutionContext.global).resource.allocated.unsafeRunSync()._1

  /*
  we should use these tests to ensure that our client class correctly handles responses from the provider - i.e. decoding, error mapping, validation
   */
  it should "get DrsHub ok status" in {
    new DrsHubClientImpl[IO](client, Uri.unsafeFromString(mockServer.getUrl))
      .fetchSystemStatus()
      .attempt
      .unsafeRunSync() shouldBe Right(true)
  }

  it should "resolve drs object" in {
    new DrsHubClientImpl[IO](client, Uri.unsafeFromString(mockServer.getUrl))
      .resolveDrsObject("drs://drs.example.com/1234567890", requestFields)
      .attempt
      .unsafeRunSync() shouldBe Right(drsResourceResponsePlaceholder)
  }
}
