package centaur
import java.util.concurrent.Executors
import cats.effect.IO
import scala.concurrent.ExecutionContext

package object TestContext {
  // This is hacky but I don't think centaur is overly concerned about exploding thread pools
  val ec = ExecutionContext.fromExecutorService(Executors.newCachedThreadPool())
  implicit val cs = IO.contextShift(ec)
}
