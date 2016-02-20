package cromwell.engine.workflow

import akka.actor.ActorSystem
import akka.testkit.{TestActorRef, ImplicitSender, TestKit}
import cromwell.backend.config.BackendConfigurationEntry
import cromwell.engine.workflow.ValidateActor.GetExecutableBackends
import cromwell.util.SampleWdl
import org.scalatest.mock.MockitoSugar
import org.scalatest.{BeforeAndAfterAll, MustMatchers, WordSpecLike}

class ValidateActorSpec extends TestKit(ActorSystem("ValidateActorSpec"))
  with WordSpecLike
  with MustMatchers
  with BeforeAndAfterAll
  with ImplicitSender
  with MockitoSugar {

  val wdlSrc = SampleWdl.HelloWorld.wdlSource()
  val wdlJson = SampleWdl.HelloWorld.wdlJson

  "ValidateActor" should {
    "successfully return a list of validated backends" in {
      val validateActor = TestActorRef(new ValidateActor(wdlSrc, Option(wdlJson), None))
      validateActor ! GetExecutableBackends
      expectMsgClass(classOf[List[BackendConfigurationEntry]])
      system.stop(validateActor)
    }
    "should fail to return a list for bad wdl" in {
      import scala.concurrent.duration._
      val validateActor = TestActorRef(new ValidateActor("Yololo!", None, None))
      validateActor ! GetExecutableBackends
      receiveWhile(1 second){
        case akka.actor.Status.Failure(something) => println(s"Message: $something")
        case oops => fail(s"This is unexpected. Msg: $oops")
      }
    }
  }

  override def afterAll(): Unit = system.shutdown()
}
