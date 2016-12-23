package cromwell.services.metadata.impl

import cromwell.services.metadata.MetadataValue
import org.scalatest.{FlatSpec, Matchers}

class MetadataValueSpec extends FlatSpec with Matchers {
  "MetadataValue" should "not NPE with a null value" in {
    val metadataValue = MetadataValue(null)
    metadataValue.value shouldBe ""
  }
}
