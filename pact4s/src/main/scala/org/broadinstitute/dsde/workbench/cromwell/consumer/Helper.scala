package org.broadinstitute.dsde.workbench.cromwell.consumer

import au.com.dius.pact.consumer.dsl.{DslPart, PactDslResponse, PactDslWithProvider}
import org.broadinstitute.dsde.workbench.model.WorkbenchEmail
import org.http4s.Credentials.Token
import org.http4s.{AuthScheme, Credentials}
import pact4s.algebras.PactBodyJsonEncoder
case object UnknownError extends Exception

object AuthHelper {
  def mockBearerHeader(workbenchEmail: WorkbenchEmail) = s"Bearer TokenFor$workbenchEmail"
  def mockAuthToken(workbenchEmail: WorkbenchEmail): Token =
    Credentials.Token(AuthScheme.Bearer, s"TokenFor$workbenchEmail")
}

object PactHelper {
  def buildInteraction(builder: PactDslResponse,
                       state: String,
                       uponReceiving: String,
                       method: String,
                       path: String,
                       requestHeaders: Seq[(String, String)],
                       status: Int,
                       responseHeaders: Seq[(String, String)],
                       body: DslPart
  ): PactDslResponse =
    builder
      .`given`(state)
      .uponReceiving(uponReceiving)
      .method(method)
      .path(path)
      .headers(scala.jdk.CollectionConverters.MapHasAsJava(requestHeaders.toMap).asJava)
      .willRespondWith()
      .status(status)
      .headers(scala.jdk.CollectionConverters.MapHasAsJava(responseHeaders.toMap).asJava)
      .body(body)

  def buildInteraction(builder: PactDslResponse,
                       state: String,
                       uponReceiving: String,
                       method: String,
                       path: String,
                       requestHeaders: Seq[(String, String)],
                       status: Int
  ): PactDslResponse =
    builder
      .`given`(state)
      .uponReceiving(uponReceiving)
      .method(method)
      .path(path)
      .headers(scala.jdk.CollectionConverters.MapHasAsJava(requestHeaders.toMap).asJava)
      .willRespondWith()
      .status(status)

  def buildInteraction(builder: PactDslResponse,
                       uponReceiving: String,
                       method: String,
                       path: String,
                       requestHeaders: Seq[(String, String)],
                       status: Int
  ): PactDslResponse =
    builder
      .uponReceiving(uponReceiving)
      .method(method)
      .path(path)
      .headers(scala.jdk.CollectionConverters.MapHasAsJava(requestHeaders.toMap).asJava)
      .willRespondWith()
      .status(status)

  def buildInteraction(builder: PactDslWithProvider,
                       state: String,
                       uponReceiving: String,
                       method: String,
                       path: String,
                       requestHeaders: Seq[(String, String)],
                       status: Int,
                       responseHeaders: Seq[(String, String)],
                       body: DslPart
  ): PactDslResponse =
    builder
      .`given`(state)
      .uponReceiving(uponReceiving)
      .method(method)
      .path(path)
      .headers(scala.jdk.CollectionConverters.MapHasAsJava(requestHeaders.toMap).asJava)
      .willRespondWith()
      .status(status)
      .headers(scala.jdk.CollectionConverters.MapHasAsJava(responseHeaders.toMap).asJava)
      .body(body)

  def buildInteraction(builder: PactDslWithProvider,
                       state: String,
                       stateParams: Map[String, Any],
                       uponReceiving: String,
                       method: String,
                       path: String,
                       requestHeaders: Seq[(String, String)],
                       status: Int,
                       responseHeaders: Seq[(String, String)],
                       body: DslPart
  ): PactDslResponse =
    builder
      .`given`(state, scala.jdk.CollectionConverters.MapHasAsJava(stateParams).asJava)
      .uponReceiving(uponReceiving)
      .method(method)
      .path(path)
      .headers(scala.jdk.CollectionConverters.MapHasAsJava(requestHeaders.toMap).asJava)
      .willRespondWith()
      .status(status)
      .headers(scala.jdk.CollectionConverters.MapHasAsJava(responseHeaders.toMap).asJava)
      .body(body)

  def buildInteraction(builder: PactDslResponse,
                       state: String,
                       stateParams: Map[String, Any],
                       uponReceiving: String,
                       method: String,
                       path: String,
                       requestHeaders: Seq[(String, String)],
                       status: Int,
                       responseHeaders: Seq[(String, String)],
                       body: DslPart
  ): PactDslResponse =
    builder
      .`given`(state, scala.jdk.CollectionConverters.MapHasAsJava(stateParams).asJava)
      .uponReceiving(uponReceiving)
      .method(method)
      .path(path)
      .headers(scala.jdk.CollectionConverters.MapHasAsJava(requestHeaders.toMap).asJava)
      .willRespondWith()
      .status(status)
      .headers(scala.jdk.CollectionConverters.MapHasAsJava(responseHeaders.toMap).asJava)
      .body(body)

  def buildInteraction(builder: PactDslResponse,
                       state: String,
                       stateParams: Map[String, Any],
                       uponReceiving: String,
                       method: String,
                       path: String,
                       requestHeaders: Seq[(String, String)],
                       status: Int,
                       responseHeaders: Seq[(String, String)],
                      ): PactDslResponse =
    builder
      .`given`(state, scala.jdk.CollectionConverters.MapHasAsJava(stateParams).asJava)
      .uponReceiving(uponReceiving)
      .method(method)
      .path(path)
      .headers(scala.jdk.CollectionConverters.MapHasAsJava(requestHeaders.toMap).asJava)
      .willRespondWith()
      .status(status)
      .headers(scala.jdk.CollectionConverters.MapHasAsJava(responseHeaders.toMap).asJava)


  def buildInteraction[A](builder: PactDslWithProvider,
                          state: String,
                          stateParams: Map[String, Any],
                          uponReceiving: String,
                          method: String,
                          path: String,
                          requestHeaders: Seq[(String, String)],
                          status: Int,
                          responseHeaders: Seq[(String, String)],
                          body: A
  )(implicit ev: PactBodyJsonEncoder[A]): PactDslResponse =
    builder
      .`given`(state, scala.jdk.CollectionConverters.MapHasAsJava(stateParams).asJava)
      .uponReceiving(uponReceiving)
      .method(method)
      .path(path)
      .headers(scala.jdk.CollectionConverters.MapHasAsJava(requestHeaders.toMap).asJava)
      .willRespondWith()
      .status(status)
      .headers(scala.jdk.CollectionConverters.MapHasAsJava(responseHeaders.toMap).asJava)
      .body(ev.toJsonString(body))

  def buildInteraction[A](builder: PactDslResponse,
                          state: String,
                          uponReceiving: String,
                          method: String,
                          path: String,
                          requestHeaders: Seq[(String, String)],
                          requestBody: A,
                          status: Int
  )(implicit ev: PactBodyJsonEncoder[A]): PactDslResponse =
    builder
      .`given`(state)
      .uponReceiving(uponReceiving)
      .method(method)
      .path(path)
      .headers(scala.jdk.CollectionConverters.MapHasAsJava(requestHeaders.toMap).asJava)
      .body(ev.toJsonString(requestBody))
      .willRespondWith()
      .status(status)

  def buildInteraction(builder: PactDslWithProvider,
                       state: String,
                       stateParams: Map[String, Any],
                       uponReceiving: String,
                       method: String,
                       path: String,
                       requestHeaders: Seq[(String, String)],
                       status: Int,
                       responseHeaders: Seq[(String, String)],
                      ): PactDslResponse =
    builder
      .`given`(state, scala.jdk.CollectionConverters.MapHasAsJava(stateParams).asJava)
      .uponReceiving(uponReceiving)
      .method(method)
      .path(path)
      .headers(scala.jdk.CollectionConverters.MapHasAsJava(requestHeaders.toMap).asJava)
      .willRespondWith()
      .status(status)
      .headers(scala.jdk.CollectionConverters.MapHasAsJava(responseHeaders.toMap).asJava)


  def buildInteraction[A](builder: PactDslResponse,
                          state: String,
                          stateParams: Map[String, Any],
                          uponReceiving: String,
                          method: String,
                          path: String,
                          requestHeaders: Seq[(String, String)],
                          requestBody: A,
                          status: Int
  )(implicit ev: PactBodyJsonEncoder[A]): PactDslResponse =
    builder
      .`given`(state, scala.jdk.CollectionConverters.MapHasAsJava(stateParams).asJava)
      .uponReceiving(uponReceiving)
      .method(method)
      .path(path)
      .headers(scala.jdk.CollectionConverters.MapHasAsJava(requestHeaders.toMap).asJava)
      .body(ev.toJsonString(requestBody))
      .willRespondWith()
      .status(status)
}
