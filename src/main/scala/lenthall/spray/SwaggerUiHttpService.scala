package lenthall.spray

import com.typesafe.config.Config
import spray.http.StatusCodes
import spray.routing.HttpService

trait SwaggerUiHttpService extends HttpService {

  val swaggerUiInfo: SwaggerUiInfo

  final def swaggerUiRoutes = {
    swaggerUiInfo match {
      case SwaggerUiInfo(swaggerUiVersion, baseUrl, docsPath, uiPath) =>
        get {
          pathPrefix(uiPath.split("/").map(segmentStringToPathMatcher).reduceLeft(_ / _)) {
            // when the user hits the doc url, redirect to the index.html with api docs specified on the url
            pathEndOrSingleSlash { context =>
              context.redirect(s"$baseUrl/$uiPath/index.html?url=$baseUrl/$docsPath", StatusCodes.TemporaryRedirect)
            } ~ getFromResourceDirectory(s"META-INF/resources/webjars/swagger-ui/$swaggerUiVersion")
          }
        }
    }
  }
}

/**
 * Attributes to render an instance of Swagger UI.
 *
 * @param swaggerUiVersion Version of the swagger-ui bundle containing the resources, for example "2.1.1".
 * @param baseUrl Base URL where the application is hosted.
 * @param docsPath Endpoint for serving the api docs. Should NOT include the baseUrl.
 * @param uiPath Path to serve the swagger-ui. Should NOT include the baseUrl.
 */
case class SwaggerUiInfo(swaggerUiVersion: String,
                         baseUrl: String = "",
                         docsPath: String = "api-docs",
                         uiPath: String = "swagger")

trait ConfigSwaggerUiHttpService extends SwaggerUiHttpService {
  def swaggerUiConfig: Config

  import lenthall.util.ScalaConfig._

  val swaggerUiInfo = {
    var info = SwaggerUiInfo(swaggerUiConfig.getString("uiVersion"))
    info = swaggerUiConfig.getStringOption("baseUrl").map(x => info.copy(baseUrl = x)) getOrElse info
    info = swaggerUiConfig.getStringOption("docsPath").map(x => info.copy(docsPath = x)) getOrElse info
    info = swaggerUiConfig.getStringOption("uiPath").map(x => info.copy(uiPath = x)) getOrElse info
    info
  }
}
