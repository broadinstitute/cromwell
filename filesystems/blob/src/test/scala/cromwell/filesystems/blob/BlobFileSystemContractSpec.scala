package cromwell.filesystems.blob

import au.com.dius.pact.consumer.PactTestExecutionContext
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

import pact4s.scalatest.RequestResponsePactForger
import au.com.dius.pact.core.model.RequestResponsePact
import au.com.dius.pact.consumer.ConsumerPactBuilder
import io.circe.Json

class BlobFileSystemContractSpec extends AnyFlatSpec with Matchers with RequestResponsePactForger {

  /**
   * we can define the folder that the pact contracts get written to upon completion of this test suite.
   */
  override val pactTestExecutionContext: PactTestExecutionContext = new PactTestExecutionContext(
    "./filesystems/blob/target/pacts"
  )

  override def pact: RequestResponsePact = ConsumerPactBuilder
      .consumer("cromwell-blob-filesystem-consumer")
      .hasPactWith("wsm-provider")
      .`given`(
        "resource exists",
        Map("id" -> resourceId, "value" -> 123) // we can use parameters to specify details about the provider state
      )
      .`given`(
        "workspace exists",
        Map("id" -> workspaceId, "value" -> 123) // we can use parameters to specify details about the provider state
      )
      .uponReceiving("Request to fetch SAS Token")
      .method("POST")
      .path(s"/api/workspaces/v1/${workspaceId}/resources/controlled/azure/storageContainer/${resourceId}/getSasToken")
      .headers("Authorization" -> AzureCredentials.getAccessToken())
      .willRespondWith()
      .status(200)
      .body(
        Json.obj("id" -> testID.asJson, "value" -> 123.asJson)
      )

}
