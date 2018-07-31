package common.util

import org.scalatest.{FlatSpec, Matchers}
import common.util.IntrospectableLazy._


class IntrospectableLazySpec extends FlatSpec with Matchers {

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

    assert(myLazy() == 4)
    assert(myLazy() == 4)
    assert(myLazy() == 4)

    assert(lazyInstantiations == 1)
    assert(myLazy.exists)
  }

}
