package cromwell.core.io

import java.nio.file.Paths

import common.assertion.CromwellTimeoutSpec
import cromwell.core.io.DefaultIoCommand.DefaultIoCopyCommand
import cromwell.core.path.DefaultPathBuilder
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers

class IoAckSpec extends AnyFlatSpecLike with CromwellTimeoutSpec with Matchers {

  "IoFailAck pattern matching" should "work for both IoFailure and IoReadForbiddenFailure" in {
    import DefaultPathBuilder._
    val command = DefaultIoCopyCommand(build(Paths.get("foo")), build(Paths.get("bar")), overwrite = false)

    val ioFailure = IoFailure(command, new RuntimeException("blah"))
    val ioReadForbiddenFailure = IoFailure(command, new RuntimeException("blah"))

    ioFailure match {
      case IoFailAck(_: IoCommand[_], _: Throwable) => ()
      case _ => fail("IoFailure did not unapply")
    }

    ioReadForbiddenFailure match {
      case IoFailAck(_: IoCommand[_], _: Throwable) => ()
      case _ => fail("IoReadForbiddenFailure did not unapply")
    }
  }
}
