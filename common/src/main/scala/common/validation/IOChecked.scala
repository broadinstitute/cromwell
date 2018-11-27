package common.validation

import cats.arrow.FunctionK
import cats.data.EitherT.fromEither
import cats.data.{EitherT, NonEmptyList, ValidatedNel}
import cats.effect.{ContextShift, IO, LiftIO}
import cats.{Applicative, Monad, Parallel}
import common.Checked
import common.exception.{AggregatedException, AggregatedMessageException, CompositeException, MessageAggregation}
import common.validation.ErrorOr.ErrorOr
import common.validation.Validation._

import scala.concurrent.Future
import scala.util.{Failure, Success, Try}

object IOChecked {
  /**
    * Type alias for EitherT[IO, NonEmptyList[String], A]
    * Represents an IO[ Either[NonEmptyList[String], A] ]
    * The monad transformer allows to flatMap over the value while keeping the IO effect as well as the list of potential errors
    */
  type IOChecked[A] = EitherT[IO, NonEmptyList[String], A]
  /**
    * Fixes the left type of Either to Throwable
    * This is useful when calling on IO[A].attempt, transforming it to IO[ Either[Throwable, A] ] == IO[ Attempt[A] ]
    */
  type Attempt[A] = Either[Throwable, A]

  /**
    * Type on which we can build an Applicative and therefore use in a Parallel[IOChecked, IOCheckedPar] (see below)
    * It uses the IO.Par type which can evaluates IOs in parallel provided a ContextShift
    */
  type IOCheckedPar[A] = IO.Par[Attempt[A]]

  /**
    * An applicative instance on Attempt which is aware of MessageAggregation and can unpack them to maintain the list
    * of error messages when combining Left values
    */
  implicit val eitherThrowableApplicative = new Applicative[Attempt] {
    override def pure[A](x: A) = Right(x)
    override def ap[A, B](ff: Attempt[A => B])(fa: Attempt[A]): Attempt[B] = {
      (fa, ff) match {
          // Both have a list or error messages, combine them in a single one
        case (Left(t1: MessageAggregation), Left(t2: MessageAggregation)) => Left(AggregatedMessageException("", t1.errorMessages ++ t2.errorMessages))
          // Only one of them is a MessageAggregation, combined the errors and the other exception in a CompositeException
        case (Left(t1: MessageAggregation), Left(t2)) => Left(CompositeException("", List(t2), t1.errorMessages))
        case (Left(t2), Left(t1: MessageAggregation)) => Left(CompositeException("", List(t2), t1.errorMessages))
          // None of them is a MessageAggregation, combine the 2 throwables in an AggregatedException
        case (Left(t1), Left(t2)) => Left(AggregatedException("", List(t1, t2)))
          // Success case, apply f on v
        case (Right(v), Right(f)) => Right(f(v))
          // Default failure case, just keep the failure
        case (Left(t1), _) => Left(t1)
        case (_, Left(t1)) => Left(t1)
      }
    }
  }

  /**
    * Derived type class of Parallel for IOChecked, using IOCheckedPar as the applicative
    * This allows to evaluate IOCheckeds in parallel by evaluating their underlying IO in Parallel.
    * We can't simply use IO.Par to evaluate them because IO.Par stops the evaluation as soon as an IO fails (actually cancels the other one)
    * This is not what we want because we want to keep the "ErrorOr" behavior where we accumulate all the errors.
    * Instead, we use IOCheckedPar which works on an Attempt[A]. The attempt corresponds to an IO[A].attempt, which will always yield a successful IO
    * but instead materialize the value as an Either[Throwable, A] (== Attempt[A]). We can then use IO.Par to evaluate all those IOs in parallel,
    * and they will all succeed. We will then need to look into the resulting Either to decide other or not the operation was actually successful or not.
    * 
    * We use MessageAggregation exceptions to pack Nel[String] into exceptions and unpack them, as we convert between IOChecked and IOCheckedPar
    * 
    * e.g:
    * import IOChecked._
    * import cats.syntax.parallel._
    * import cats.instances.list._
    * 
    * val ec = ExecutionContext.fromExecutor(Executors.newFixedThreadPool(5))
    * implicit val contextShift = IO.contextShift(ec)
    * 
    * // All the IOs will evaluate in parallel on the above execution context
    * List(0, 1, 2).parTraverse[IOChecked, IOCheckedPar, Char](i => 
    *   EitherT.liftF(IO(Right(i.toChar))) // == IOChecked[Char]
    * )
    */
  implicit def ioCheckedParallel(implicit cs: ContextShift[IO]): Parallel[IOChecked, IOCheckedPar] = new Parallel[IOChecked, IOCheckedPar] {
    // Applicative instance for IOCheckedPar is the one for Io.Par composed with Attempt, since IOCheckedPar[A] == IO.Par[Attempt[A]]
    override def applicative: Applicative[IOCheckedPar] = IO.parApplicative(cs).compose[Attempt]
    // Monad instance for IOChecked can be summoned from cats implicit (IOChecked is simply an EitherT which already has a cats defined Monad instance)
    override def monad: Monad[IOChecked] = implicitly[Monad[IOChecked]]

    /**
      * How to convert an IOCheckedPar to an IOChecked
      * To do this we first unpack the IO.Par to an IO using IO.Par.unwrap
      * Then we examine the result of the IO by flatMapping on it
      * If the value is a Left (a failure), and the failure is a MessageAggregation, we need to extract those failure messages
      *   into an Nel[String]
      * Otherwise it's simply a failed or successful IO
      * We now have an IO[ Either[NonEmptyList[String], A] ] which we can wrap into an IOChecked with EitherT.apply
      */
    override def sequential: FunctionK[IOCheckedPar, IOChecked] =
      new FunctionK[IOCheckedPar, IOChecked] { 
        def apply[A](fa: IOCheckedPar[A]): IOChecked[A] = EitherT {
          IO.Par.unwrap(fa: IO.Par[Attempt[A]]) flatMap {
            case Left(t: MessageAggregation) => IO.pure(Left(NonEmptyList.fromListUnsafe(t.errorMessages.toList)))
            case Left(t) => IO.raiseError(t)
            case Right(v) => IO.pure(Right(v))
          }
        }
      }

    /**
      * How to convert an IOChecked to an IOCheckedPar
      * Here we get the IO value of the IOChecked and call attempt on it
      * We go from an IO[ Either[NonEmptyList[String], A] ] to an IO[ Either[Throwable, Either[NonEmptyList[String], A] ] ]
      * If we convert that IO to an IO.Par using IO.Par.apply, we get a 
      * IO.Par[ Either[Throwable, Either[NonEmptyList[String], A] ] ] 
      *   == IO.Par[Attempt[ Either[NonEmptyList[String], A] ] ]
      *   == IOCheckedPar[ Either[NonEmptyList[String], A] ]
      * To get an IOCheckedPar[A], we convert the 
      *  Either[Throwable, Either[NonEmptyList[String], A] ] to Either[Throwable, A ]
      *  by creating a AggregatedMessageException from the Nel[String]
      * This gives us an IO.Par[ Either[Throwable, A] ] == IO.Par[ Attempt[A] ] == IOCheckedPar[A]
      */
    override def parallel: FunctionK[IOChecked, IOCheckedPar] =
      new FunctionK[IOChecked, IOCheckedPar] {
        override def apply[A](fa: IOChecked[A]): IOCheckedPar[A] = IO.Par {
          fa.value.attempt map {
            case Left(t) => Left(t)
            case Right(Left(t)) => Left(AggregatedMessageException("", t.toList))
            case Right(Right(value)) => Right(value)
          }
        }
      }
  }

  def error[A](error: String, tail: String*): IOChecked[A] = EitherT.leftT {
    NonEmptyList.of(error, tail: _*)
  }
  
  def pure[A](value: A): IOChecked[A] = EitherT.pure(value)

  def goIOChecked[A](f: => A): IOChecked[A] = Try(f).toIOChecked

  implicit val IOCheckedLiftIO = new LiftIO[IOChecked] {
    override def liftIO[A](ioa: IO[A]): IOChecked[A] = EitherT.liftF[IO, NonEmptyList[String], A](ioa)
  }

  implicit class EnhancedIOChecked[A](val p: IOChecked[A]) extends AnyVal {
    import cats.syntax.either._

    def toChecked: Checked[A] = {
      Try(p.value.unsafeRunSync()) match {
        case Success(r) => r
        case Failure(f) => NonEmptyList.one(f.getMessage).asLeft
      }
    }

    def toErrorOr: ErrorOr[A] = toChecked.toValidated

    def unsafeToEither(): Either[NonEmptyList[String], A] = p.value.unsafeRunSync()

    def unsafe(context: String) = unsafeToEither().unsafe(context)

    def contextualizeErrors(context: String): IOChecked[A] = {
      p.leftMap({ errors =>
        val total = errors.size
        errors.zipWithIndex map { case (e, i) => s"Failed to $context (reason ${i + 1} of $total): $e" }
      })
    }
  }

  implicit class TryIOChecked[A](val t: Try[A]) extends AnyVal {
    def toIOChecked: IOChecked[A] = t.toErrorOr.toIOChecked
  }
  
  implicit class FutureIOChecked[A](val future: Future[A]) extends AnyVal {
    def toIOChecked: IOChecked[A] = IO.fromFuture(IO(future)).to[IOChecked]
  }

  implicit class ErrorOrIOChecked[A](val e: ErrorOr[A]) extends AnyVal {
    def toIOChecked: IOChecked[A] = fromEither[IO](e.toEither)
  }

  implicit class CheckedIOChecked[A](val c: Checked[A]) extends AnyVal {
    def toIOChecked: IOChecked[A] = fromEither[IO](c)
  }

  implicit class ValidIOChecked[A](val obj: A) extends AnyVal {
    def validIOChecked: IOChecked[A] = IOChecked.pure(obj)
  }

  implicit class InvalidIOChecked(val obj: String) extends AnyVal {
    def invalidIOChecked[A]: IOChecked[A] = error(obj)
  }

  implicit class OptionIOChecked[A](val o: Option[A]) extends AnyVal {
    def toIOChecked(errorMessage: String): IOChecked[A] = {
      EitherT.fromOption(o, NonEmptyList.of(errorMessage))
    }
  }

  type IOCheckedValidated[A] = IO[ValidatedNel[String, A]]
}
