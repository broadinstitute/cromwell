package cromwell.filesystems.blob

import au.com.dius.pact.consumer.PactTestExecutionContext
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

import pact4s.scalatest.RequestResponsePactForger
import au.com.dius.pact.core.model.RequestResponsePact

class BlobFileSystemContractSpec extends AnyFlatSpec with Matchers with RequestResponsePactForger {

  override def pact: RequestResponsePact = ???

  /*
  we can define the folder that the pact contracts get written to upon completion of this test suite.
   */
  override val pactTestExecutionContext: PactTestExecutionContext = new PactTestExecutionContext(
    "./example/resources/pacts"
  )
}
