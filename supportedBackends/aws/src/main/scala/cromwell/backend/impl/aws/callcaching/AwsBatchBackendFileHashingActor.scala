/*
 * Copyright 2018 Amazon.com, Inc. or its affiliates.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *
 *  1. Redistributions of source code must retain the above copyright notice,
 *  this list of conditions and the following disclaimer.
 *
 *  2. Redistributions in binary form must reproduce the above copyright
 *  notice, this list of conditions and the following disclaimer in the
 *  documentation and/or other materials provided with the distribution.
 *
 *  3. Neither the name of the copyright holder nor the names of its
 *  contributors may be used to endorse or promote products derived from
 *  this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
 *  BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
 *  FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL
 *  THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 *  INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 *  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 *  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 *  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 *  STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
 *  IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */
package cromwell.backend.impl.aws.callcaching

import cromwell.backend.standard.callcaching.{StandardFileHashingActor, StandardFileHashingActorParams}
import cromwell.filesystems.s3.batch.S3BatchCommandBuilder
import cromwell.backend.BackendInitializationData
import cromwell.backend.impl.aws.AwsBatchBackendInitializationData
import cromwell.backend.impl.aws.AWSBatchStorageSystems
import cromwell.core.callcaching.FileHashStrategy
import cromwell.core.io.DefaultIoCommandBuilder

class AwsBatchBackendFileHashingActor(standardParams: StandardFileHashingActorParams)
    extends StandardFileHashingActor(standardParams) {

  override val defaultHashingStrategies: Map[String, FileHashStrategy] = Map(
    ("s3", FileHashStrategy.ETag)
  )

  override val ioCommandBuilder = BackendInitializationData
    .as[AwsBatchBackendInitializationData](standardParams.backendInitializationDataOption)
    .configuration
    .batchAttributes
    .fileSystem match {
    case AWSBatchStorageSystems.s3 => S3BatchCommandBuilder
    case _ => DefaultIoCommandBuilder
  }
}
