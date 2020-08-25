package common.util

import common.util.IntrospectableLazy._
import org.scalatest.concurrent.Futures
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

import scala.concurrent.ExecutionContext.Implicits.global
import scala.concurrent.duration._
import scala.concurrent.{Await, Future}

class IntrospectableLazySpec extends AnyFlatSpec with Matchers with Futures {

  behavior of "IntrospectableLazy"

  it should "work the way it's designed" in {
    var lazyInstantiations = 0

    def lazyContents = {
      lazyInstantiations += 1
      4
    }

    val myLazy = lazily { lazyContents }

    assert(lazyInstantiations == 0)
    assert(!myLazy.exists)

    // Fails without `synchronized { ... }`
    Await.result(Future.sequence(
      Seq.fill(100)(Future {
        myLazy() shouldBe 4
      })
    ), 1.seconds)

    assert(lazyInstantiations == 1)
    assert(myLazy.exists)
  }

}
