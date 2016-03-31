package lenthall.exception

import java.io.{PrintWriter, PrintStream}

import org.scalatest.{FlatSpecLike, Matchers}

class ExceptionAggregationSpec extends FlatSpecLike with Matchers{

  "MessageAggregation" should "aggregate messages" in {
    val aggregatedException = new Exception with MessageAggregation {
      override def exceptionContext: String = "Bouuhhh"
      override def errorMessages: Traversable[String] = List("Didn't work", "didn't work either")
    }

    aggregatedException.getMessage shouldBe
      """Bouuhhh
        |Didn't work
        |didn't work either""".stripMargin
  }

  "ThrowableAggregation" should "aggregate throwables" in {
    val fakeStackTrace1 = Array(new StackTraceElement("SomeClass", "Some method", "Some file", 5))
    val fakeStackTrace2 = Array(new StackTraceElement("SomeOtherClass", "Some other method", "Some other file", 12))

    val exception1 = new RuntimeException("What is wrong with you ?")
    exception1.setStackTrace(fakeStackTrace1)
    val exception2 = new IllegalArgumentException("Please stop coding altogether")
    exception2.setStackTrace(fakeStackTrace1)

    val throwableAggregation = new Exception with ThrowableAggregation {
      override def exceptionContext: String = "Clearly not working"
      override def throwables: Traversable[Throwable] = List(exception1, exception2)
    }

    val aggregatedException = new AggregatedException("Clearly not working", List(exception1, exception2))

    List(throwableAggregation, aggregatedException) foreach { e =>
      e.getMessage shouldBe
        """Clearly not working
          |What is wrong with you ?
          |Please stop coding altogether""".stripMargin

      e.getSuppressed should contain theSameElementsAs List(exception1, exception2)
    }
  }
}
