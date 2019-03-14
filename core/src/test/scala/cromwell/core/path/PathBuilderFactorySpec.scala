package cromwell.core.path

import akka.actor.ActorSystem
import com.typesafe.config.{Config, ConfigFactory}
import cromwell.core.{TestKitSuite, WorkflowOptions}
import org.scalatest.concurrent.ScalaFutures
import org.scalatest.{FlatSpecLike, Matchers}

import scala.concurrent.{ExecutionContext, Future}

class PathBuilderFactorySpec extends TestKitSuite with FlatSpecLike with ScalaFutures with Matchers {
  behavior of "PathBuilderFactory"
  implicit val ec = system.dispatcher
  
  it should "sort factories when instantiating path builders" in {
    val factory1 = new MockPathBuilderFactory(ConfigFactory.empty(), ConfigFactory.parseString("name=factory1"))
    val factory2 = new MockPathBuilderFactory(ConfigFactory.empty(), ConfigFactory.parseString("name=factory2"))
    PathBuilderFactory
      .instantiatePathBuilders(List(DefaultPathBuilderFactory, factory1, factory2), WorkflowOptions.empty).map({ pathBuilders =>
      pathBuilders.last shouldBe DefaultPathBuilder
      // check that the order of the other factories has not been changed
      pathBuilders.map(_.name) shouldBe List("factory1", "factory2", DefaultPathBuilder.name)
    }).futureValue
  }
}

class MockPathBuilderFactory(globalConfig: Config, val instanceConfig: Config) extends cromwell.core.path.PathBuilderFactory {
  override def withOptions(options: WorkflowOptions)(implicit as: ActorSystem, ec: ExecutionContext) = Future.successful(
    new PathBuilder {
      override def name = instanceConfig.getString("name")
      override def build(pathAsString: String) = throw new UnsupportedOperationException
    }
  )
}
