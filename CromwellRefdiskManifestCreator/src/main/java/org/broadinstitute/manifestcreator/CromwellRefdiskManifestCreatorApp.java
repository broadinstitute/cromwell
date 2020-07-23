package org.broadinstitute.manifestcreator;

import com.fasterxml.jackson.databind.ObjectMapper;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;
import org.broadinstitute.manifestcreator.exception.CRC32Exception;
import org.broadinstitute.manifestcreator.model.ReferenceDiskManifest;
import org.broadinstitute.manifestcreator.model.ReferenceFile;

import java.io.*;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Stack;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.stream.Stream;
import java.util.zip.CRC32;

public class CromwellRefdiskManifestCreatorApp {

  private static final Logger logger = LogManager.getLogger(CromwellRefdiskManifestCreatorApp.class);

  public static void main(String[] args) throws IOException, InterruptedException {
    Configurator.setRootLevel(Level.DEBUG);

    if (args.length < 4) {
      logger.error("Wrong number of parameters.");
      logger.error("Usage: java -jar cromwell-refdisk-manifest-creator-app.jar <number of parallel threads> " +
              "<image_name> <directory_path_to_scan> <output file path>");
      System.exit(1);
    }

    String directoryToScan = args[2];
    File rootDirectory = new File(directoryToScan);
    if (!rootDirectory.exists() || !rootDirectory.isDirectory()) {
      logger.error("Root directory " + directoryToScan + " doesn't exist.");
      System.exit(1);
    }

    int totalNumberOfFiles = countFiles(rootDirectory);
    logger.info("Parsing " + totalNumberOfFiles + " files");

    ReferenceDiskManifest manifest = new ReferenceDiskManifest();
    String imageName = args[1];
    manifest.setImageIdentifier(imageName);

    int nThreads = Integer.parseInt(args[0]);
    CountDownLatch countDownLatch = new CountDownLatch(nThreads);
    ExecutorService executorService = Executors.newFixedThreadPool(nThreads);

    Path rootPath = rootDirectory.toPath();
    Stack<File> stack = new Stack<>();
    stack.push(rootDirectory);
    List<ReferenceFile> referenceFiles = Collections.synchronizedList(new ArrayList<>());
    for (int i = 0; i < nThreads; i++) {
      executorService.submit(() -> {
        try {
          doWork(rootPath, referenceFiles, stack, totalNumberOfFiles, countDownLatch);
        } catch (CRC32Exception e) {
          logger.error("Error occurred. Shutting down.", e);
          executorService.shutdownNow();
          System.exit(1);
        }
      });
    }

    countDownLatch.await();

    String manifestFilePath = args[3];
    File manifestFile = new File(manifestFilePath);
    if (manifestFile.exists() || !manifestFile.createNewFile()) {
      logger.error("File " + manifestFile.getAbsolutePath() + " already exists or cannot be created.");
    } else {
      referenceFiles.forEach(referenceFile -> manifest.getFiles().add(referenceFile));
      new ObjectMapper().writerWithDefaultPrettyPrinter().writeValue(manifestFile, manifest);
      logger.info("Finished. Total number of reference files: " + manifest.getFiles().size());
    }

    executorService.shutdown();
  }

  private static void doWork(Path rootPath, List<ReferenceFile> accumulator, Stack<File> stack, int totalNumberOfFiles, CountDownLatch countDownLatch) throws CRC32Exception {
    while(!stack.isEmpty() || accumulator.size() != totalNumberOfFiles) {
      if (stack.isEmpty()) {
        continue;
      }

      File curFile = stack.pop();

      if (curFile.isDirectory()) {
        File[] children = curFile.listFiles();
        if (children != null) {
          Stream.of(children).forEach(stack::push);
        }
      } else {
        ReferenceFile refFile = new ReferenceFile();

        String relativePath = rootPath.relativize(curFile.toPath()).toString();
        refFile.setPath(relativePath);

        long crc32 = calculateCrc32(curFile);
        refFile.setCrc32(crc32);

        accumulator.add(refFile);
        logger.info(Thread.currentThread().getName() + " finished processing file " + curFile.getAbsolutePath());
      }
    }

    logger.info("Thread " + Thread.currentThread().getName() + " finished processing.");
    countDownLatch.countDown();
  }

  private static long calculateCrc32(File file) throws CRC32Exception {
    CRC32 crc32 = new CRC32();
    int bufferSize = 32768;
    try (InputStream is = new BufferedInputStream(new FileInputStream(file), bufferSize)) {
      byte[] buff = new byte[bufferSize];
      int bytesRead;
      while ((bytesRead = is.read(buff)) != -1) {
        crc32.update(buff, 0, bytesRead);
      }
    } catch (IOException e) {
      throw new CRC32Exception("Cannot read from file " + file.getAbsolutePath(), e);
    }
    return crc32.getValue();
  }

  private static int countFiles(File directory) {
    int count = 0;
    Stack<File> stack = new Stack<>();
    stack.push(directory);
    while(!stack.isEmpty()) {
      File curFile = stack.pop();
      if(curFile.isDirectory()) {
        File[] children = curFile.listFiles();
        if (children != null) {
          Stream.of(children).forEach(stack::push);
        }
      } else {
        count++;
      }
    }
    return count;
  }

}
