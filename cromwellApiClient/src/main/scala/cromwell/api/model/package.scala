package cromwell.api

import java.time.OffsetDateTime

import akka.http.scaladsl.model.HttpResponse
import akka.http.scaladsl.unmarshalling.Unmarshal
import akka.stream.ActorMaterializer
import cats.arrow.FunctionK
import cats.data.EitherT
import cats.effect.{ContextShift, IO, Timer}
import cromwell.api.CromwellClient.UnsuccessfulRequestException
import cromwell.api.model.TimeUtil._
import spray.json.{DefaultJsonProtocol, JsString, JsValue, RootJsonFormat}

import scala.concurrent.duration.FiniteDuration
import scala.concurrent.{ExecutionContext, Future}

package object model {

  implicit val OffsetDateTimeJsonFormat = OffsetDateTimeJsonFormatter.OffsetDateTimeFormat

  object OffsetDateTimeJsonFormatter extends DefaultJsonProtocol {
    object OffsetDateTimeFormat extends RootJsonFormat[OffsetDateTime] {
      def write(offsetDateTime: OffsetDateTime) = new JsString(offsetDateTime.toUtcMilliString)
      def read(value: JsValue) = value match {
        case JsString(string) => OffsetDateTime.parse(string)
        case other => throw new UnsupportedOperationException(s"Cannot deserialize $other into an OffsetDateTime")
      }
    }
  }

  type FailureResponseOrT[A] = EitherT[IO, HttpResponse, A]
  val FailureResponseOrT = EitherT

  implicit class EnhancedFutureHttpResponse(val responseFuture: Future[HttpResponse]) extends AnyVal {
    def asFailureResponseOrT: FailureResponseOrT[HttpResponse] = {
      val ioHttpResponse = IO.fromFuture(IO(responseFuture))
      val ioEither = ioHttpResponse map {
        case response if response.status.isFailure() => Left(response)
        case response => Right(response)
      }
      FailureResponseOrT(ioEither)
    }
  }

  implicit class EnhancedFailureResponseOrHttpResponseT(val responseIoT: FailureResponseOrT[HttpResponse])
    extends AnyVal {
    def asHttpResponse: Future[HttpResponse] = {
      val io = responseIoT.value map {
        case Left(response) => response
        case Right(response) => response
      }
      io.unsafeToFuture()
    }
  }

  implicit class EnhancedFailureResponseOrT[SuccessType](val responseIoT: FailureResponseOrT[SuccessType]) extends AnyVal {
    final def timeout(duration: FiniteDuration)
                     (implicit timer: Timer[IO], cs: ContextShift[IO]): FailureResponseOrT[SuccessType] = {
      EitherT(responseIoT.value.timeout(duration))
    }

    def asIo(implicit materializer: ActorMaterializer, executionContext: ExecutionContext): IO[SuccessType] = {
      responseIoT.value flatMap {
        case Left(response) =>
          IO.fromFuture(IO {
            Unmarshal(response.entity).to[String] flatMap { responseString =>
              Future.failed(UnsuccessfulRequestException(responseString, response))
            }
          })
        case Right(a) => IO.pure(a)
      }
    }

    /**
      * Transforms the IO error from one type to another.
      */
    def mapErrorWith(mapper: Throwable => IO[Nothing]): FailureResponseOrT[SuccessType] = {
      def handleErrorIo[A](ioIn: IO[A]): IO[A] = {
        ioIn handleErrorWith mapper
      }

      responseIoT.mapK(FunctionK.lift(handleErrorIo))
    }
  }
}
