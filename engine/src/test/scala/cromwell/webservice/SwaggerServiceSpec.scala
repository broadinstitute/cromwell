package cromwell.webservice

import akka.http.scaladsl.model.StatusCodes
import akka.http.scaladsl.testkit.ScalatestRouteTest
import io.swagger.models.properties.RefProperty
import io.swagger.parser.SwaggerParser
import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.{FlatSpec, Matchers}
import org.yaml.snakeyaml.constructor.Constructor
import org.yaml.snakeyaml.error.YAMLException
import org.yaml.snakeyaml.nodes.MappingNode
import org.yaml.snakeyaml.{Yaml => SnakeYaml}

import scala.collection.JavaConverters._

class SwaggerServiceSpec extends FlatSpec with SwaggerService with ScalatestRouteTest with Matchers
  with TableDrivenPropertyChecks {
  def actorRefFactory = system

  behavior of "SwaggerService"

  it should "return valid yaml" in {
    Get("/swagger/cromwell.yaml") ~>
      swaggerUiResourceRoute ~>
      check {
        status should be(StatusCodes.OK)

        val body = responseAs[String]
        val yaml = new SnakeYaml(new UniqueKeyConstructor()).loadAs(body, classOf[java.util.Map[String, AnyRef]])

        yaml.get("swagger") should be("2.0")
      }
  }

  it should "return valid swagger" in {
    Get("/swagger/cromwell.yaml") ~>
      swaggerUiResourceRoute ~>
      check {
        status should be(StatusCodes.OK)

        // https://github.com/swagger-api/swagger-parser/issues/976
        val swaggerBugMsg = "needs to be defined as a path parameter in path or operation level"

        val body = responseAs[String]
        val resultWithInfo = new SwaggerParser().readWithInfo(body)
        val swaggerVersion = resultWithInfo.getSwagger.getSwagger
        val swaggerMessages = resultWithInfo.getMessages.asScala.filterNot(_ contains swaggerBugMsg)

        swaggerVersion should be("2.0")
        swaggerMessages should be(empty)

        resultWithInfo.getSwagger.getDefinitions.asScala foreach {
          // If no properties, `getProperties` returns `null` instead of an empty map
          case (defKey, defVal) => Option(defVal.getProperties).map(_.asScala).getOrElse(Map.empty) foreach {
            /*
            Two against one.
            Swagger parser implementation lets a RefProperty have descriptions.
            http://swagger.io/specification/#referenceObject & http://editor.swagger.io both say it's ref ONLY!
             */
            case (propKey, propVal: RefProperty) =>
              withClue(s"RefProperty $defKey.$propKey has a description: ") {
                propVal.getDescription should be(null)
              }
            case _ => /* ignore */
          }
        }
      }
  }

  it should "return the index.html" in {
    Get("/swagger/index.html?url=/swagger/cromwell.yaml") ~>
      swaggerUiResourceRoute ~>
      check {
        assertResult(StatusCodes.OK) {
          status
        }
        assertResult("<!-- HTML for s") {
          responseAs[String].take(15)
        }
      }
  }

  it should "return status OK when getting OPTIONS on paths" in {
    val pathExamples = Table("path", "/", "/swagger", "/swagger/cromwell.yaml", "/swagger/index.html", "/api",
      "/api/workflows/", "/api/workflows/v1", "/workflows/v1/outputs", "/workflows/v1/status",
      "/api/workflows/v1/validate", "/workflows", "/workflows/v1", "/workflows/v1/outputs", "/workflows/v1/status",
      "/workflows/v1/validate")

    forAll(pathExamples) { path =>
      Options(path) ~>
        swaggerUiResourceRoute ~>
        check {
          assertResult(StatusCodes.OK) {
            status
          }
          assertResult("OK") {
            responseAs[String]
          }
        }
    }
  }
}

/**
  * SnakeYaml only support detecting duplicate keys after jumping over the hurdles.
  *
  * https://bitbucket.org/asomov/snakeyaml/issues/337
  * https://groups.google.com/forum/#!topic/snakeyaml-core/M-TRKg-0-F4
  * https://code.google.com/archive/p/snakeyaml/issues/139
  *
  * Adapted from:
  * https://bitbucket.org/asomov/snakeyaml/src/e9cd9f5e8d76c61eb983e29b3dc039c1fac9c393/src/test/java/org/yaml/snakeyaml/issues/issue139/UniqueKeyTest.java?fileviewer=file-view-default#UniqueKeyTest.java-43:62
  */
class UniqueKeyConstructor extends Constructor {

  import java.util.{Map => JMap}

  override protected def constructMapping2ndStep(node: MappingNode, mapping: JMap[AnyRef, AnyRef]): Unit = {
    val nodeValue = node.getValue
    for (tuple <- nodeValue.asScala) {
      val keyNode = tuple.getKeyNode
      val valueNode = tuple.getValueNode
      val key = constructObject(keyNode)
      if (key != null) {
        key.hashCode(); // check circular dependencies
      }
      val value = constructObject(valueNode)
      val old = mapping.put(key, value)
      if (old != null) {
        throw new YAMLException(s"The key is not unique $key:\n[old value]\n$old\n[new value]\n$value")
      }
    }
  }
}
