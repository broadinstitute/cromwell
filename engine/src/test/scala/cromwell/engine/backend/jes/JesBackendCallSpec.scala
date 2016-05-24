package cromwell.engine.backend.jes

import akka.actor.ActorSystem
import cromwell.engine.backend.{WorkflowDescriptor, BackendCallJobDescriptor}
import cromwell.engine.backend.runtimeattributes.CromwellRuntimeAttributes
import cromwell.engine.workflow.BackendCallKey
import org.scalatest.{FlatSpec, Matchers}
import org.specs2.mock.Mockito
import wdl4s.CallInputs

/** FIXME this class references a JesBackendCall which is soon to disappear, move all this to JesBackendSpec. */
class JesBackendCallSpec extends FlatSpec with Matchers with Mockito {
  "JesBackendCall" should "return preemptible = true only in the correct cases" in {

    val backendCallKeyWithAttempt1 = mock[BackendCallKey]
    backendCallKeyWithAttempt1.attempt returns 1

    val backendCallKeyWithAttempt2 = mock[BackendCallKey]
    backendCallKeyWithAttempt2.attempt returns 2

    val workflow = mock[WorkflowDescriptor]
    val backend = new JesBackend(ActorSystem("Jessie"))
    workflow.backend returns backend

    class MaxMockingDescriptor(max: Int, key: BackendCallKey) extends BackendCallJobDescriptor(workflow, key, mock[CallInputs]) {
      val attributes = mock[CromwellRuntimeAttributes]
      override lazy val callRuntimeAttributes: CromwellRuntimeAttributes = attributes.preemptible returns max
    }

    class Attempt1(max: Int) extends MaxMockingDescriptor(max, backendCallKeyWithAttempt1)
    class Attempt2(max: Int) extends MaxMockingDescriptor(max, backendCallKeyWithAttempt2)

    val descriptorWithMax0AndKey1 = new Attempt1(max = 0)
    descriptorWithMax0AndKey1.preemptible shouldBe false

    val descriptorWithMax1AndKey1 = new Attempt1(max = 1)
    descriptorWithMax1AndKey1.preemptible shouldBe true

    val descriptorWithMax2AndKey1 = new Attempt1(max = 2)
    descriptorWithMax2AndKey1.preemptible shouldBe true

    val descriptorWithMax1AndKey2 = new Attempt2(max = 1)
    descriptorWithMax1AndKey2.preemptible shouldBe false

    val descriptorWithMax2AndKey2 = new Attempt2(max = 2)
    descriptorWithMax2AndKey2.preemptible shouldBe true
  }
}
