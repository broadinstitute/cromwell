package cromwell.services.instrumentation

import cats.data.NonEmptyList
import cromwell.core.TestKitSuite
import cromwell.services.instrumentation.CromwellInstrumentation.InstrumentationPath
import org.scalatest.PrivateMethodTester
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers

class CromwellInstrumentationSpec extends TestKitSuite with AnyFlatSpecLike with Matchers with PrivateMethodTester {

  it should "keep order from constructor" in {
    InstrumentationPath
      .withParts("a", "b", "c", "d").internalPath shouldBe NonEmptyList
      .of(Left("a"), Left("b"), Left("c"), Left("d"))
  }

  it should "allow path to start with a high-variant part" in {
    InstrumentationPath
      .withHighVariantPart("label" -> "a")
      .withParts("b", "c").internalPath shouldBe NonEmptyList
      .of(Right("label" -> "a"), Left("b"), Left("c"))
  }

  private val pathOne = InstrumentationPath
    .withParts("a")
    .withHighVariantPart("label-b" -> "b")
    .withParts("c", "d")
    .withHighVariantPart("label-e", "e")
    .withParts(List("f", "g"))

  it should "keep order across different parts" in {
    pathOne.internalPath shouldBe NonEmptyList
      .of(Left("a"), Right("label-b" -> "b"), Left("c"), Left("d"), Right("label-e" -> "e"), Left("f"), Left("g"))
  }

  it should "provide an 'all parts in order' view like the old Nel representation" in {
    pathOne.getPath shouldBe NonEmptyList
      .of("a", "b", "c", "d", "e", "f", "g")
  }

  it should "provide a 'low-variant parts in order with others in map' view" in {
    pathOne.getPathLowVariants shouldBe ((
      NonEmptyList.of("a", "c", "d", "f", "g"),
      Map("label-b" -> "b", "label-e" -> "e")
    ))
  }

  private val pathTwo = InstrumentationPath
    .withHighVariantPart("label-a" -> "A")
    .withHighVariantPart("label-b" -> "B")
    .withHighVariantPart("label-c" -> "C")

  it should "handle no normal parts" in {
    pathTwo.getPath shouldBe NonEmptyList
      .of("A", "B", "C")
  }

  it should "use a high-variant part if the low-variant parts would be empty" in {
    pathTwo.getPathLowVariants shouldBe ((
      NonEmptyList.of("A"),
      Map("label-b" -> "B", "label-c" -> "C")
    ))
  }

  it should "handle concatenation of paths" in {
    pathOne.concat(pathTwo).internalPath shouldBe NonEmptyList
      .of(Left("a"), Right("label-b" -> "b"), Left("c"), Left("d"), Right("label-e" -> "e"), Left("f"), Left("g"),
        Right("label-a" -> "A"), Right("label-b" -> "B"), Right("label-c" -> "C"))
  }

  it should "provide 'all parts in order' view across concatenations" in {
    pathOne.concat(pathTwo).getPath shouldBe NonEmptyList
      .of("a", "b", "c", "d", "e", "f", "g", "A", "B", "C")
  }

  it should "make labels last-distinct when providing 'low-variant' view" in {
    pathOne.concat(pathTwo).getPathLowVariants shouldBe ((
      NonEmptyList.of("a", "c", "d", "f", "g"),
      Map("label-b" -> "B", "label-e" -> "e", "label-a" -> "A", "label-c" -> "C")
    ))
  }

  it should "handle status codes" in {
    InstrumentationPath
      .withParts("something")
      .withStatusCodeFailure(Some(200))
      .internalPath shouldBe NonEmptyList
      .of(Left("something"), Right("code" -> "200"))
    InstrumentationPath
      .withParts("something")
      .withStatusCodeFailure(None)
      .internalPath shouldBe NonEmptyList
      .of(Left("something"))
  }

  it should "handle throwables" in {
    InstrumentationPath
      .withParts("something")
      .withThrowable(new RuntimeException(), _ => Some(404))
      .internalPath shouldBe NonEmptyList
      .of(Left("something"), Right("code" -> "404"))
    InstrumentationPath
      .withParts("something")
      .withThrowable(new RuntimeException(), _ => None)
      .internalPath shouldBe NonEmptyList
      .of(Left("something"))
  }
}
