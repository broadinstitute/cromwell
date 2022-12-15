package cromwell.filesystems.blob

import au.com.dius.pact.consumer.{ConsumerPactBuilder, PactTestExecutionContext}
import au.com.dius.pact.core.model.RequestResponsePact
import io.circe.Json
import io.circe.syntax._
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import pact4s.circe.implicits._
import pact4s.scalatest.RequestResponsePactForger

class BlobFileSystemContractSpec extends AnyFlatSpec with Matchers with RequestResponsePactForger {

  val resourceId = "";
  val workspaceId = "";
  /**
   * we can define the folder that the pact contracts get written to upon completion of this test suite.
   */
  override val pactTestExecutionContext: PactTestExecutionContext = new PactTestExecutionContext(
    "./filesystems/blob/target/pacts"
  )

  /**
   * The following is an outline for what a RequestResponsePact might look like for requesting
   * from WSM form the blob filesystem. To generate a pact file however tests must be written
   * against this pact to validate that the agreement is met. These tests get run with other
   * scala tests on build and if the tests pass when run a pact file will be generated locally
   */
  override def pact: RequestResponsePact = ConsumerPactBuilder
      .consumer("cromwell-blob-filesystem-consumer")
      .hasPactWith("wsm-provider")
      .`given`(
        "resource exists",
        Map("id" -> resourceId.asJson, "value" -> 123.asJson) // we can use parameters to specify details about the provider state
      )
      .`given`(
        "workspace exists",
        Map("id" -> workspaceId, "value" -> 123) // we can use parameters to specify details about the provider state
      )
      .uponReceiving("Request to fetch SAS Token")
      .method("POST")
      .path(s"/api/workspaces/v1/${workspaceId}/resources/controlled/azure/storageContainer/${resourceId}/getSasToken")
      .headers("Authorization" -> "sampleToken")
      .willRespondWith()
      .status(200)
      .body(
        Json.obj("id" -> "".asJson, "value" -> 123.asJson)
      ).toPact()

}
