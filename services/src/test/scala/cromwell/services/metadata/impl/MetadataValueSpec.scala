package cromwell.services.metadata.impl

import cromwell.services.metadata.MetadataValue
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

class MetadataValueSpec extends AnyFlatSpec with Matchers {
  "MetadataValue" should "not NPE with a null value" in {
    val metadataValue = MetadataValue(null)
    metadataValue.value shouldBe ""
  }
}
