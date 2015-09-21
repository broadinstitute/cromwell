package cromwell.webservice

import com.gettyimages.spray.swagger.SwaggerHttpService
import com.typesafe.config.Config
import com.wordnik.swagger.model.ApiInfo
import spray.http.StatusCodes
import spray.routing.HttpService

// TODO: This swagger boilerplate, and the Http.Bind in CromwellServer, should be in a dsde-common
trait SwaggerConfigHttpService extends SwaggerHttpService with SwaggerUiHttpService {
  def swaggerConfig: Config

  override lazy val baseUrl: String = swaggerConfig.getString("baseUrl")
  // NOTE: These two values in the conf must be kept in sync with changes to build.sbt
  // Validated by SwaggerConfigHttpServiceSpec
  override lazy val apiVersion: String = swaggerConfig.getString("api.version")
  override lazy val swaggerUiVersion: String = swaggerConfig.getString("ui.version")

  override lazy val docsPath: String = optString("docsPath", super.docsPath)
  override lazy val swaggerVersion: String = optString("swaggerVersion", super.swaggerVersion)
  override lazy val uiPath: String = optString("ui.path", super.uiPath)

  override lazy val apiInfo: Option[ApiInfo] =
    if (swaggerConfig.hasPath("api.info"))
      Option(new ApiInfo(
        swaggerConfig.getString("api.info"),
        swaggerConfig.getString("api.description"),
        swaggerConfig.getString("api.termsOfServiceUrl"),
        swaggerConfig.getString("api.contact"),
        swaggerConfig.getString("api.license"),
        swaggerConfig.getString("api.licenseUrl")
      ))
    else
      None

  private def optString(path: String, default: => String) = {
    if (swaggerConfig.hasPath(path))
      swaggerConfig.getString(path)
    else
      default
  }
}

trait SwaggerUiHttpService extends HttpService {
  this: SwaggerHttpService =>

  /**
   * @return path to serve the swagger-ui
   */
  def uiPath = "swagger"

  /**
   * The swagger-ui version of the swagger-ui bundle containing the resources, for example "2.1.8-M1"
   *
   * @return swagger-ui version
   */
  def swaggerUiVersion: String

  final def uiRoutes = {
    get {
      pathPrefix(uiPath) {
        // when the user hits the doc url, redirect to the index.html with api docs specified on the url
        pathEndOrSingleSlash { context =>
          context.redirect(s"$baseUrl/$uiPath/index.html?url=$baseUrl/swagger/cromwell.yaml", StatusCodes.TemporaryRedirect)
        } ~ {
          getFromResourceDirectory("swagger/") ~ getFromResourceDirectory(s"META-INF/resources/webjars/swagger-ui/$swaggerUiVersion")
        }
      }
    }
  }
}
