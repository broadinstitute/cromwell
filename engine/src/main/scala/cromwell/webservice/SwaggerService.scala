package cromwell.webservice

import cromwell.webservice.routes.CromwellApiService

trait SwaggerService extends SwaggerUiResourceHttpService {
  override def swaggerServiceName = "cromwell"

  override def swaggerUiVersion = CromwellApiService.swaggerUiVersion
}
