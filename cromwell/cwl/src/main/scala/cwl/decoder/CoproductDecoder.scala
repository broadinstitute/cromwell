package cwl.decoder

import io.circe._
import shapeless.{:+:, CNil, Coproduct, Inl, Inr}
import cats.syntax.show._

/**
  * Custom Coproduct decoder that keeps the last error message if all of the decoding attempts fail.
  */
trait CoproductDecoder {

  implicit final val decodeCNil: Decoder[CNil] = new Decoder[CNil] {
    def apply(c: HCursor): Decoder.Result[CNil] = Left(DecodingFailure("CNil", c.history))
  }

  implicit final val encodeCNil: Encoder[CNil] = new Encoder[CNil] {
    def apply(a: CNil): Json = sys.error("Cannot encode CNil")
  }

  implicit final def decodeCCons[L, R <: Coproduct](implicit
                                                    decodeL: Decoder[L],
                                                    decodeR: Decoder[R]
                                                   ): Decoder[L :+: R] = {
    decodeL.map[L :+: R](Inl(_)).handleErrorWith(error => decodeR.withErrorMessage(error.show).map(Inr(_)))
  }

  implicit final def encodeCCons[L, R <: Coproduct](implicit
                                                    encodeL: Encoder[L],
                                                    encodeR: Encoder[R]
                                                   ): Encoder[L :+: R] = new Encoder[L :+: R] {
    def apply(a: L :+: R): Json = a match {
      case Inl(l) => encodeL(l)
      case Inr(r) => encodeR(r)
    }
  }
}
