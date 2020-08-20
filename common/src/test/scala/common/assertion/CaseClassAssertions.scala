package common.assertion

import org.scalatest.matchers.should.Matchers

object CaseClassAssertions extends Matchers {
  implicit class ComparableCaseClass[A <: Product](actualA: A) {
    // Assumes that expectedA and actualA are the same type. If we don't subtype case classes, that should hold...
    def shouldEqualFieldwise(expectedA: A): Unit = {
      (0 until actualA.productArity) foreach { i =>
        actualA.productElement(i) should be(expectedA.productElement(i))
      }
    }
  }
}
