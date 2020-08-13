package org.broadinstitute.manifestcreator;

import com.fasterxml.jackson.databind.ObjectMapper;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;
import org.broadinstitute.manifestcreator.exception.CRC32CException;
import org.broadinstitute.manifestcreator.model.ReferenceDiskManifest;
import org.broadinstitute.manifestcreator.model.ReferenceFile;

import java.io.*;
import java.nio.file.Path;
import java.util.*;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.stream.Stream;
import java.util.zip.CRC32C;

import static java.util.Collections.emptyList;
import static java.util.stream.Collectors.partitioningBy;

public class CromwellRefdiskManifestCreatorApp {

  private static final Logger logger = LogManager.getLogger(CromwellRefdiskManifestCreatorApp.class);

  public static void main(String[] args) throws IOException, InterruptedException {
    Configurator.setRootLevel(Level.INFO);

    Arguments inputArguments = parseArguments(args);

    ReferenceDiskManifest manifest = createManifestForDirectory(inputArguments);

    logger.info("Writing file to disk...");
    File manifestFile = new File(inputArguments.manifestFilePath);
    if (manifestFile.exists() || !manifestFile.createNewFile()) {
      logger.error("File {} already exists or cannot be created.", manifestFile.getAbsolutePath());
      printUsageAndExit();
    } else {
      new ObjectMapper().writerWithDefaultPrettyPrinter().writeValue(manifestFile, manifest);
      logger.info("Completed");
    }
  }

  static ReferenceDiskManifest createManifestForDirectory(Arguments inputArguments) throws InterruptedException {
    File rootDirectory = new File(inputArguments.directoryToScan);
    if (!rootDirectory.exists() || !rootDirectory.isDirectory()) {
      logger.error("Root directory {} doesn't exist.", inputArguments.directoryToScan);
      printUsageAndExit();
    }

    logger.info("Populating list of files to process...");
    Stack<File> allFilesInRootRecurStack = prepopulateStackWithFiles(rootDirectory);
    logger.info("Need to process {} files", allFilesInRootRecurStack.size());

    ReferenceDiskManifest manifest = new ReferenceDiskManifest();
    manifest.setImageIdentifier(inputArguments.imageName);
    manifest.setDiskSizeGb(inputArguments.diskSizeGb);

    CountDownLatch countDownLatch = new CountDownLatch(inputArguments.nThreads);
    ExecutorService executorService = Executors.newFixedThreadPool(inputArguments.nThreads);

    Path rootPath = rootDirectory.toPath();
    List<ReferenceFile> referenceFiles = Collections.synchronizedList(new ArrayList<>());
    for (int i = 0; i < inputArguments.nThreads; i++) {
      executorService.submit(() -> {
        try {
          processFileStackIntoReferences(rootPath, referenceFiles, allFilesInRootRecurStack, countDownLatch);
        } catch (CRC32CException e) {
          logger.error("Error occurred. Shutting down.", e);
          System.exit(1);
        }
      });
    }

    countDownLatch.await();
    executorService.shutdown();

    manifest.getFiles().addAll(referenceFiles);
    logger.info("Finished created manifest. Total number of reference files: {}", manifest.getFiles().size());

    return manifest;
  }

  private static void processFileStackIntoReferences(Path rootPath,
                             List<ReferenceFile> accumulator,
                             Stack<File> fileStack,
                             CountDownLatch countDownLatch) throws CRC32CException {
    List<ReferenceFile> interimResult = new ArrayList<>();
    while(!fileStack.isEmpty()) {
      File curFile;
      try {
        curFile = fileStack.pop();
      } catch (EmptyStackException e) {
        // it's ok - benign race condition
        break;
      }

      String relativePath = rootPath.relativize(curFile.toPath()).toString();
      ReferenceFile refFile = new ReferenceFile();
      refFile.setPath(relativePath);

      long crc32c = calculateCrc32c(curFile);
      refFile.setCrc32c(crc32c);

      interimResult.add(refFile);
      logger.info( "{} finished processing file {}", Thread.currentThread().getName(), curFile.getAbsolutePath());
    }
    accumulator.addAll(interimResult);

    logger.info("Thread {} finished processing.", Thread.currentThread().getName());
    countDownLatch.countDown();
  }

  private static long calculateCrc32c(File file) throws CRC32CException {
    CRC32C crc32c = new CRC32C();
    int bufferSize = 32768;
    try (InputStream is = new BufferedInputStream(new FileInputStream(file), bufferSize)) {
      byte[] buff = new byte[bufferSize];
      int bytesRead;
      while ((bytesRead = is.read(buff)) != -1) {
        crc32c.update(buff, 0, bytesRead);
      }
    } catch (IOException e) {
      throw new CRC32CException("Cannot read from file " + file.getAbsolutePath(), e);
    }
    return crc32c.getValue();
  }

  private static Stack<File> prepopulateStackWithFiles(File directory) {
    Stack<File> dirStack = new Stack<>();
    Stack<File> fileStack = new Stack<>();
    dirStack.push(directory);
    while(!dirStack.isEmpty()) {
      File curFile = dirStack.pop();
      if(curFile.isDirectory()) {
        File[] children = curFile.listFiles();
        if (children != null) {
          Map<Boolean, List<File>> dirsAndFiles = Stream.of(children).collect(partitioningBy(File::isDirectory));
          dirStack.addAll(dirsAndFiles.getOrDefault(true, emptyList()));
          fileStack.addAll(dirsAndFiles.getOrDefault(false, emptyList()));
        }
      } else {
        fileStack.add(curFile);
      }
    }
    return fileStack;
  }

  private static Arguments parseArguments(String[] args) {
    if (args.length < 5) {
      logger.error("Wrong number of parameters.");
      printUsageAndExit();
    }

    int nThreads = parseIntArg(args[0], "number of threads");
    String imageName = args[1];
    int diskSizeGb = parseIntArg(args[2], "disk size");
    String directoryToScan = args[3];
    String manifestFilePath = args[4];

    return new Arguments(nThreads, imageName, diskSizeGb, directoryToScan, manifestFilePath);
  }

  private static int parseIntArg(String arg, String argName) {
    try {
      int intVal = Integer.parseInt(arg);
      if (intVal < 1) {
        throw new NumberFormatException();
      }
      return intVal;
    } catch (NumberFormatException e) {
      logger.error("Wrong argument type - {} should be a positive integer", argName);
      printUsageAndExit();
      return -1;
    }
  }

  private static void printUsageAndExit() {
    logger.error("Usage: java -jar cromwell-refdisk-manifest-creator-app.jar <number of parallel threads> " +
            "<image identifier> <disk size Gb> <directory path to scan> <output file path>");
    System.exit(1);
  }

  static class Arguments {

    int nThreads;
    String imageName;
    int diskSizeGb;
    String directoryToScan;
    String manifestFilePath;

    Arguments(int nThreads, String imageName, int diskSizeGb, String directoryToScan, String manifestFilePath) {
      this.nThreads = nThreads;
      this.imageName = imageName;
      this.diskSizeGb = diskSizeGb;
      this.directoryToScan = directoryToScan;
      this.manifestFilePath = manifestFilePath;
    }

  }
}
