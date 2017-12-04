package cromwell.webservice

trait SwaggerService extends SwaggerUiResourceHttpService {
  override def swaggerServiceName = "cromwell"

  override def swaggerUiVersion = CromwellApiService.swaggerUiVersion
}
