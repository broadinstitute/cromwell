package cromwell.filesystems.demo.dos

import com.fasterxml.jackson.databind.node.ArrayNode
import com.fasterxml.jackson.databind.{JsonNode, ObjectMapper}
import com.typesafe.config.Config
import com.typesafe.scalalogging.StrictLogging
import net.ceedubs.ficus.Ficus._
import net.thisptr.jackson.jq.{JsonQuery, Scope}
import org.apache.http.client.methods.HttpPost
import org.apache.http.entity.{ContentType, StringEntity}
import org.apache.http.impl.client.{CloseableHttpClient, HttpClientBuilder}
import org.apache.http.util.EntityUtils
import org.apache.http.{HttpResponse, HttpStatus}

import scala.collection.JavaConverters._
import scala.util.Try

case class DemoDosResolver(config: Config) extends StrictLogging {

  private val GcsScheme = "gs"
  private val DosPathToken = s"$${dosPath}"
  private lazy val objectMapper = new ObjectMapper()
  private lazy val marthaUri = config.getString("drs.martha.url")
  private lazy val marthaRequestJsonTemplate = config.getString("drs.martha.request.json-template")
  private lazy val marthaResponseJqFilterString = config.getString("drs.martha.response.jq-filter")
  private lazy val marthaResponseJqJsonQuery = JsonQuery.compile(marthaResponseJqFilterString)
  private lazy val marthaResponseJqScope = {
    val scope = Scope.newEmptyScope
    scope.loadFunctions(Thread.currentThread().getContextClassLoader)
    scope
  }
  private lazy val marthaDebug = config.getOrElse("demo.dos.martha.debug", false)

  private def debugResponse(response: HttpResponse): String = if (marthaDebug) s"\n$response" else ""

  protected def createClient(): CloseableHttpClient = HttpClientBuilder.create().build()

  def getContainerRelativePath(demoDosPath: DemoDosPath): String = {
    val client = createClient()
    try {
      val post = new HttpPost(marthaUri)
      val requestJson = marthaRequestJsonTemplate.replace(DosPathToken, demoDosPath.pathAsString)
      post.setEntity(new StringEntity(requestJson, ContentType.APPLICATION_JSON))

      if (marthaDebug) {
        logger.info(s"Martha Request:\n$post\n${post.getEntity}\n${EntityUtils.toString(post.getEntity)}")
      }

      val response = client.execute(post)
      val responseEntityOption = Option(response.getEntity).map(EntityUtils.toString)

      if (marthaDebug) {
        logger.info(s"Martha Response:\n$response${responseEntityOption.mkString("\n", "", "")}")
      }

      def throwUnexpectedResponse: Nothing = {
        throw new RuntimeException(
          s"Unexpected response looking up ${demoDosPath.pathAsString} from $marthaUri.${debugResponse(response)}")
      }

      try {
        val content = if (response.getStatusLine.getStatusCode == HttpStatus.SC_OK) {
          responseEntityOption.getOrElse(throwUnexpectedResponse)
        } else {
          throwUnexpectedResponse
        }

        val contentNode = objectMapper.readTree(content)
        val filteredNodes = marthaResponseJqJsonQuery(marthaResponseJqScope, contentNode).asScala
        val resolvedUrlOption = filteredNodes.toStream.map(getResolvedUrl).collectFirst {
          case Some(url) => url
        }
        val pathWithoutSchemeOption = resolvedUrlOption.map(_.substring(GcsScheme.length + 3))

        pathWithoutSchemeOption getOrElse throwUnexpectedResponse
      } finally {
        Try(response.close())
        ()
      }
    } finally {
      Try(client.close())
      ()
    }
  }

  private def getResolvedUrl(node: JsonNode): Option[String] = {
    node match {
      case _ if node.isMissingNode => None
      case _ if node.isArray =>
        val array = node.asInstanceOf[ArrayNode]
        array.elements.asScala.map(getResolvedUrl) collectFirst {
          case Some(element) => element
        }
      case _ if node.isTextual =>
        Option(node.asText) flatMap { text =>
          if (text.startsWith(GcsScheme + "://")) Option(text) else None
        }
      case _ => None
    }
  }
}
