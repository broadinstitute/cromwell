package lenthall.legacy

import org.scalatest.{FlatSpec, Matchers}
import lenthall.legacy.TwoElevenSupport._

import scala.util.{Failure, Success}

class TacticallyEnhancedTryTest extends FlatSpec with Matchers {
  behavior of "TacticallyEnhancedTry"

  it should "convert 'Success'es right" in {
    Success("hello").tacticalToEither should be(Right("hello"))
  }

  it should "'Failure's should just be left" in {
    val error = new Exception("*that* was the bad thing")
    Failure(error).tacticalToEither should be(Left(error))
  }
}
