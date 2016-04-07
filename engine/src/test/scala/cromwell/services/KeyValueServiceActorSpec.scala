package cromwell.services

import java.util.UUID

import akka.actor.ActorSystem
import akka.pattern.ask
import akka.util.Timeout
import com.typesafe.config.ConfigFactory
import cromwell.CromwellTestkitSpec.TestWorkflowManagerSystem
import cromwell.engine.backend.local.LocalBackend
import cromwell.engine.db.slick.SlickDataAccess
import cromwell.engine.workflow.BackendCallKey
import cromwell.engine.{WorkflowDescriptor, WorkflowId, WorkflowSourceFiles}
import cromwell.services.KeyValueServiceActor.{Get, KeyValuePair, Put}
import org.scalatest.time.{Seconds, Millis, Span}
import org.scalatest.{FlatSpec, Matchers}
import wdl4s.Call
import scala.concurrent.duration._

import scala.concurrent.Await

class KeyValueServiceActorSpec extends FlatSpec with Matchers {
  val workflowManagerSystem = new TestWorkflowManagerSystem
  val actorSystem = workflowManagerSystem.actorSystem
  val databaseConfig = ConfigFactory.parseString(
    s"""
       |db.url = "jdbc:hsqldb:mem:$${slick.uniqueSchema};shutdown=false;hsqldb.tx=mvcc"
       |db.driver = "org.hsqldb.jdbcDriver"
       |driver = "slick.driver.HsqldbDriver$$"
       |schema.manager = slick
       |""".stripMargin)
  val dataAccess = new SlickDataAccess(databaseConfig)
  val workflowId = WorkflowId(UUID.randomUUID())
  val localBackend = LocalBackend(workflowManagerSystem.actorSystem)

  val sources = WorkflowSourceFiles(
    wdlSource="""task a {command{}}
                |workflow w {
                |  call a
                |  call a as b
                |  call a as c
                |}
              """.stripMargin,
    inputsJson="{}",
    workflowOptionsJson="{}"
  )

  val descriptor = WorkflowDescriptor(id = workflowId, sourceFiles = sources)
  def call(name: String): Call = descriptor.namespace.workflow.calls.find(_.unqualifiedName == name).get
  val callA = call("a")

  implicit val timeout = Timeout(5.seconds)
  implicit val ec = actorSystem.dispatcher

  it should "foobar" in {
    val kv = actorSystem.actorOf(
      KeyValueServiceActor.props(dataAccess, workflowId, "", ConfigFactory.parseString(""))
    )

    val bc = BackendCallKey(callA, None, 1)

    val x = for {
      _ <- dataAccess.createWorkflow(descriptor, Nil, descriptor.namespace.workflow.calls, localBackend)
      p <- ask(kv, Put(bc, "key", "value"))
      _ = println("foobar")
      g <- ask(kv, Get(bc, "key")).mapTo[KeyValuePair]
    } yield g

    val r = Await.result(x, 2.seconds)
    println(r)
    r shouldEqual KeyValuePair("key", Some("value"))
  }
}
