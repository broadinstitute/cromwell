package cromwell.engine.workflow.tokens

import cromwell.core.JobExecutionToken.JobExecutionTokenType
import cromwell.engine.workflow.tokens.UnhoggableTokenPool.{HogLimitExceeded, TokenHoggingLease, TokenTypeExhausted, TokensAvailable}
import org.scalatest.concurrent.Eventually
import org.scalatest.{FlatSpec, Matchers}

import scala.concurrent.duration._

class UnhoggableTokenPoolSpec extends FlatSpec with Matchers with Eventually {

  override val patienceConfig = PatienceConfig(timeout = scaled(5.seconds), interval = scaled(1.second))
  implicit val patience = patienceConfig

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
      val pool = new UnhoggableTokenPool(tokenType)
      pool.hogLimitOption should be(expectedHogLimit)
    }
  }

  hogLimitingTokenTypeToHogLimit foreach { case (tokenType, hogLimit) =>
    it should s"enforce the hogLimit of $hogLimit for ${tokenType.hogFactor} groups for $tokenType" in {
      val pool = new UnhoggableTokenPool(tokenType)

      (0 until tokenType.hogFactor) foreach { hogGroupNumber =>
        val hogGroup = s"hogGroup$hogGroupNumber"
        // Step one: hog as many tokens as the group is able!:
        (0 until hogLimit) foreach { index =>
          pool.tryAcquire(hogGroup) match {
            case _: TokenHoggingLease => // great!
            case TokenTypeExhausted => fail(s"Unhoggable token pool ran out after $index tokens distributed to $hogGroupNumber")
            case HogLimitExceeded => fail(s"Unhoggable token pool making unfounded accusations of hogging after $index tokens distributed to $hogGroupNumber")
          }

          val acquiredTokensForGroup = index + 1

          val poolState = pool.poolState
          poolState.leased should be(hogGroupNumber * hogLimit + acquiredTokensForGroup)

          val hogGroupState = poolState.hogGroups.get.find(_.hogGroup == hogGroup).get
          hogGroupState.used should be(acquiredTokensForGroup)
          hogGroupState.atLimit should be(acquiredTokensForGroup == hogLimit)
        }

        // Step two: one more is too many...
        val nextAcquisition = pool.tryAcquire(hogGroup)
        val poolState = pool.poolState

        if (hogGroupNumber != tokenType.hogFactor - 1) {
          // While the overall token type is not full we get the hog limit push-back:
          nextAcquisition should be(HogLimitExceeded)
          pool.available("other") should be(UnhoggableTokenPool.TokensAvailable)
          poolState.available should be(true)

        } else {
          // The final time, the entire token type is used up so we get this instead:
          nextAcquisition should be(TokenTypeExhausted)
          pool.available("other") should be(UnhoggableTokenPool.TokenTypeExhausted)
          poolState.available should be(false)
        }
      }
    }
  }

  it should "allow tokens to be returned" in {
    // A pool distributing tokens with a hogLimit of 2:
    val hogLimitPool = new UnhoggableTokenPool(JobExecutionTokenType("backend", Some(150), 75))

    // Use all the "group1" tokens:
    val lease1 = hogLimitPool.tryAcquire("group1").asInstanceOf[TokenHoggingLease]
    val lease2 = hogLimitPool.tryAcquire("group1").asInstanceOf[TokenHoggingLease]
    hogLimitPool.available("group1") shouldBe HogLimitExceeded
    hogLimitPool.tryAcquire("group1") should be(HogLimitExceeded)

    lease1.release()
    eventually { hogLimitPool.available("group1") shouldBe TokensAvailable }
    hogLimitPool.tryAcquire("group1") match {
      case _: TokenHoggingLease => // Great!
      case other => fail(s"expected lease but got $other")
    }
    hogLimitPool.tryAcquire("group1") should be(HogLimitExceeded)

    lease2.release()
    eventually { hogLimitPool.available("group1") shouldBe TokensAvailable }
    hogLimitPool.tryAcquire("group1") match {
      case _: TokenHoggingLease => // Great!
      case other => fail(s"expected lease but got $other")
    }
    hogLimitPool.tryAcquire("group1") should be(HogLimitExceeded)
  }

}
