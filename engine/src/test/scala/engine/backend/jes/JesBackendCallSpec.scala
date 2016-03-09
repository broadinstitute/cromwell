package cromwell.engine.backend.jes

import java.net.URL

import cromwell.engine.WorkflowDescriptor
import cromwell.engine.backend.BackendCallJobDescriptor
import cromwell.engine.workflow.BackendCallKey
import org.scalatest.{FlatSpec, Matchers}
import org.slf4j.{Logger, LoggerFactory}
import org.specs2.mock.Mockito
import wdl4s.CallInputs

class JesBackendCallSpec extends FlatSpec with Matchers with Mockito {
  "JesBackendCall" should "return preemptible = true only in the correct cases" in {
    val logger: Logger = LoggerFactory.getLogger("JesBackendCallSpec logger")
    val project = "project"
    val bucket = "gs://bucket"
    val endpoint =  new URL("http://endpoint")

    val backendCallKeyWithAttempt1 = mock[BackendCallKey]
    backendCallKeyWithAttempt1.attempt returns 1

    val backendCallKeyWithAttempt2 = mock[BackendCallKey]
    backendCallKeyWithAttempt2.attempt returns 2

    val callDescriptor1 = BackendCallJobDescriptor(mock[WorkflowDescriptor], backendCallKeyWithAttempt1, mock[CallInputs])
    val callDescriptor2 = BackendCallJobDescriptor(mock[WorkflowDescriptor], backendCallKeyWithAttempt2, mock[CallInputs])

    val backendCallWithMax0AndKey1 = new JesBackendCall(mock[JesBackend], callDescriptor1, None) {
      override lazy val maxPreemption = 0
    }
    backendCallWithMax0AndKey1.preemptible shouldBe false

    val backendCallWithMax1AndKey1 = new JesBackendCall(mock[JesBackend], callDescriptor1, None) {
      override lazy val maxPreemption = 1
    }
    backendCallWithMax1AndKey1.preemptible shouldBe true

    val backendCallWithMax2AndKey1 = new JesBackendCall(mock[JesBackend], callDescriptor1, None) {
      override lazy val maxPreemption = 2
    }
    backendCallWithMax2AndKey1.preemptible shouldBe true

    val backendCallWithMax1AndKey2 = new JesBackendCall(mock[JesBackend], callDescriptor2, None) {
      override lazy val maxPreemption = 1
    }
    backendCallWithMax1AndKey2.preemptible shouldBe false

    val backendCallWithMax2AndKey2 = new JesBackendCall(mock[JesBackend], callDescriptor2, None) {
      override lazy val maxPreemption = 2
    }
    backendCallWithMax2AndKey2.preemptible shouldBe true
  }
}
