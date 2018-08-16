package cromwell.engine.workflow.tokens

import cromwell.core.JobExecutionToken.JobExecutionTokenType
import cromwell.engine.workflow.tokens.UnhoggableTokenPool.{ComeBackLater, Oink, TokenHoggingLease}
import org.scalatest.{FlatSpec, Matchers}

class UnhoggableTokenPoolSpec extends FlatSpec with Matchers {

  val hogLimitingTokenTypeToHogLimit = List(
    JobExecutionTokenType("backend", Some(150), 2) -> 75,
    JobExecutionTokenType("backend", Some(150), 75) -> 2,
    JobExecutionTokenType("backend", Some(150), 150) -> 1
  )

  val tokenTypeToHogLimit = List(
    JobExecutionTokenType("backend", Some(150), 1) -> None,
    JobExecutionTokenType("backend", Some(150), 200) -> Some(1),
    JobExecutionTokenType("backend", None, 1) -> None,
    JobExecutionTokenType("backend", None, 150) -> None
  ) ++ (hogLimitingTokenTypeToHogLimit map { case (k,v) => (k, Some(v)) })

  tokenTypeToHogLimit foreach { case (tokenType, expectedHogLimit) =>
    it should s"correctly calculate hogLimit for $tokenType as $expectedHogLimit" in {
      val pool = UnhoggableTokenPool(tokenType)
      pool.hogLimitOption should be(expectedHogLimit)
    }
  }

  hogLimitingTokenTypeToHogLimit foreach { case (tokenType, hogLimit) =>
    it should s"enforce the hogLimit of $hogLimit for ${tokenType.hogFactor} groups for $tokenType" in {
      val pool = UnhoggableTokenPool(tokenType)

      (0 until tokenType.hogFactor) foreach { hogGroupNumber =>
        val hogGroup = s"hogGroup$hogGroupNumber"
        // Step one: hog as many tokens as the group is able!:
        (0 until hogLimit) foreach { index =>
          pool.tryAcquire(hogGroup) match {
            case _: TokenHoggingLease => // great!
            case ComeBackLater => fail(s"Unhoggable token pool ran out after $index tokens distributed to $hogGroupNumber")
            case Oink => fail(s"Unhoggable token pool making unfounded accusations of hogging after $index tokens distributed to $hogGroupNumber")
          }
        }
        // Step two: one more is too many...
        pool.tryAcquire(hogGroup) should be(Oink)
      }

      // Any more requests should be told to come back later:
      pool.tryAcquire("somebody else") should be(ComeBackLater)
    }
  }

  it should "allow tokens to be returned" in {
    val hogLimit2Pool = UnhoggableTokenPool(JobExecutionTokenType("backend", Some(150), 75))

    // Use all the "group1" tokens:
    val lease1 = hogLimit2Pool.tryAcquire("group1").asInstanceOf[TokenHoggingLease]
    val lease2 = hogLimit2Pool.tryAcquire("group1").asInstanceOf[TokenHoggingLease]
    hogLimit2Pool.tryAcquire("group1") should be(Oink)

    lease1.release()
    hogLimit2Pool.tryAcquire("group1") match {
      case _: TokenHoggingLease => // Great!
      case other => fail(s"expected lease but got $other")
    }
    hogLimit2Pool.tryAcquire("group1") should be(Oink)

    lease2.release()
    hogLimit2Pool.tryAcquire("group1") match {
      case _: TokenHoggingLease => // Great!
      case other => fail(s"expected lease but got $other")
    }
    hogLimit2Pool.tryAcquire("group1") should be(Oink)
  }

}
