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
import scala.util.Try
import cromwell.backend.standard.callcaching.StandardFileHashingActor.SingleFileHashRequest
import cromwell.core.path.DefaultPathBuilder

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
  // get backend config.
  val aws_config = BackendInitializationData.as[AwsBatchBackendInitializationData](standardParams.backendInitializationDataOption).configuration
  
  // custom strategy to handle efs (local) files, in case sibling-md5 file is present. 
  override def customHashStrategy(fileRequest: SingleFileHashRequest): Option[Try[String]] = {
    val file = DefaultPathBuilder.get(fileRequest.file.valueString)
    if (aws_config.efsMntPoint.isDefined && file.toString.startsWith(aws_config.efsMntPoint.getOrElse("--")) && aws_config.checkSiblingMd5.getOrElse(false)) {
            val md5 = file.sibling(s"${file.toString}.md5")
            // check existance of the file : 
            if (!file.exists) {
                // if missing, cache hit is invalid; return invalid md5
                Some("File Missing").map(str => Try(str))
            }
            // check existence of the sibling file
            else if (md5.exists) {
                // read the file.
                val md5_value: Option[String] = Some(md5.contentAsString.split("\\s+")(0))
                md5_value.map(str => Try(str))
            } else {
                // File present, but no sibling found, fall back to default.
                None
            }
            
    } else {
        // Detected non-EFS file: return None
        None  
    }
  }
}
