package drs.localizer.downloaders

import cats.effect.{ExitCode, IO}

trait Downloader {

  def download: IO[ExitCode]

}
