package lenthall.exception

import java.io.FileNotFoundException
import java.nio.file.NoSuchFileException

import org.scalatest.{FlatSpecLike, Matchers}

class ExceptionAggregationSpec extends FlatSpecLike with Matchers{

  "MessageAggregation" should "aggregate messages" in {
    val aggregatedException = new Exception with MessageAggregation {
      override def exceptionContext: String = "Bouuhhh"
      override def errorMessages: Traversable[String] = List("Didn't work", "didn't work either")
    }

    aggregatedException.getMessage shouldBe
      """Bouuhhh:
        |Didn't work
        |didn't work either""".stripMargin
  }

  "AggregatedMessageException" should "aggregate empty messages" in {
    val aggregatedMessageException = AggregatedMessageException("Bouuhhh", Traversable.empty)
    aggregatedMessageException.getMessage shouldBe "Bouuhhh"
  }

  "ThrowableAggregation" should "aggregate throwables" in {
    val fakeStackTrace1 = Array(new StackTraceElement("SomeClass", "Some method", "Some file", 5))

    val exception1 = new RuntimeException("What is wrong with you ?")
    exception1.setStackTrace(fakeStackTrace1)
    val exception2 = new IllegalArgumentException("Please stop coding altogether")
    exception2.setStackTrace(fakeStackTrace1)

    val throwableAggregation = new Exception with ThrowableAggregation {
      override def exceptionContext: String = "Clearly not working"
      override def throwables: Traversable[Throwable] = List(exception1, exception2)
    }

    val aggregatedException = AggregatedException("Clearly not working", List(exception1, exception2))

    List(throwableAggregation, aggregatedException) foreach { e =>
      e.getMessage shouldBe
        """Clearly not working:
          |What is wrong with you ?
          |Please stop coding altogether""".stripMargin

      e.getSuppressed should contain theSameElementsAs List(exception1, exception2)
    }
  }
  
  "ThrowableAggregation" should "aggregate throwable aggregations recursively" in {
    val exception1 = new RuntimeException("Nope")
    val exception2 = new RuntimeException("Still nope")
    val subAggregatedException = AggregatedException("Nope exception", List(exception1, exception2))
    val exception3 = new RuntimeException("Yep Exception")
    val aggregatedException = AggregatedException("This is why nothing works", List(subAggregatedException, exception3))
    
    aggregatedException.getMessage shouldBe """This is why nothing works:
                                              |Nope exception:
                                              |	Nope
                                              |	Still nope
                                              |Yep Exception""".stripMargin
  }

  "ThrowableAggregation" should "aggregate causes" in {
    val exception1 = new RuntimeException("Nope", new RuntimeException("because of nope"))
    val aggregatedException = AggregatedException("Nope exception", List(exception1))

    aggregatedException.getMessage shouldBe """Nope exception:
                                              |Nope
                                              |	because of nope""".stripMargin
  }

  "ThrowableAggregation" should "add cause to file not found exceptions" in {
    val exception1 = new FileNotFoundException("fileA")
    val exception2 = new NoSuchFileException("fileB")
    val aggregatedException = AggregatedException("Trying to read some files", List(exception1, exception2))

    aggregatedException.getMessage shouldBe """Trying to read some files:
                                              |File not found fileA
                                              |File not found fileB""".stripMargin
  }
}
