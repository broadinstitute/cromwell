package cromwell.services.metadata.impl

import common.assertion.CromwellTimeoutSpec
import cromwell.services.metadata.{MetadataBoolean, MetadataInt, MetadataNumber, MetadataString, MetadataValue}
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

class MetadataValueSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers {
  "MetadataValue" should "not NPE with a null value" in {
    val metadataValue = MetadataValue(null)
    metadataValue.value shouldBe ""
  }

  it should "have the correct data type for metadata values" in {
    val metadataInt = MetadataValue(100)
    metadataInt.valueType shouldBe MetadataInt

    val metadataFloat = MetadataValue(5.5)
    metadataFloat.valueType shouldBe MetadataNumber

    val metadataBoolean = MetadataValue(true)
    metadataBoolean.valueType shouldBe MetadataBoolean

    val metadataString = MetadataValue("hello world")
    metadataString.valueType shouldBe MetadataString
  }
}
