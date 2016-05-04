package cromwell.webservice

import io.swagger.models.properties.RefProperty
import io.swagger.parser.SwaggerParser
import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.{FlatSpec, Matchers}
import org.yaml.snakeyaml.constructor.Constructor
import org.yaml.snakeyaml.error.YAMLException
import org.yaml.snakeyaml.nodes.MappingNode
import org.yaml.snakeyaml.{Yaml => SnakeYaml}
import spray.http._
import spray.testkit.ScalatestRouteTest

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

        /*
        BUG: Right now swagger-parser says that type is unexpected in security definitions. Should be fixed in
        https://github.com/swagger-api/swagger-parser/pull/232/files#diff-392413c3c16ae4930f710f815575830dR30

        Still don't know why they choose not to print the proper location.
         */
        val swaggerBugMsg = "attribute type is unexpected"

        val body = responseAs[String]
        val resultWithInfo = new SwaggerParser().readWithInfo(body)

        resultWithInfo.getSwagger.getSwagger should be("2.0")
        resultWithInfo.getMessages.asScala.filterNot(_ == swaggerBugMsg) should be(empty)

        resultWithInfo.getSwagger.getDefinitions.asScala foreach {
          case (defKey, defVal) => defVal.getProperties.asScala foreach {
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
    Get("/swagger/index.html") ~>
      swaggerUiResourceRoute ~>
      check {
        assertResult(StatusCodes.OK) {
          status
        }
        assertResult("<!DOCTYPE html>") {
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
