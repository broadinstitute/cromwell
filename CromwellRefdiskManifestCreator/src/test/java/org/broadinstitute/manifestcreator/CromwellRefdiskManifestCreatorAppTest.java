package org.broadinstitute.manifestcreator;

import com.fasterxml.jackson.databind.ObjectMapper;
import org.broadinstitute.manifestcreator.CromwellRefdiskManifestCreatorApp.Arguments;
import org.broadinstitute.manifestcreator.model.ReferenceDiskManifest;
import org.junit.Test;

import java.io.File;
import java.io.IOException;

import static org.hamcrest.MatcherAssert.*;
import static org.hamcrest.Matchers.*;

public class CromwellRefdiskManifestCreatorAppTest {

  @Test
  public void testManifestCreation() throws InterruptedException, IOException {
    int nThreads = 1;
    String imageName = "testImageName";
    String dirToScan = "src/test/resources/test-directory";
    int diskSizeGb = 500;
    Arguments args = new Arguments(nThreads, imageName, diskSizeGb, dirToScan, "");
    ReferenceDiskManifest actualManifest = CromwellRefdiskManifestCreatorApp.createManifestForDirectory(args);

    String pathToExpectedManifest = "src/test/resources/reference_manifest.json";
    ReferenceDiskManifest expectedManifest = new ObjectMapper().readValue(new File(pathToExpectedManifest), ReferenceDiskManifest.class);

    assertThat(actualManifest.getImageIdentifier(), is(equalTo(expectedManifest.getImageIdentifier())));
    assertThat(actualManifest.getDiskSizeGb(), is(equalTo(expectedManifest.getDiskSizeGb())));
    assertThat(actualManifest.getFiles().size(), is(equalTo(expectedManifest.getFiles().size())));
    assertThat(actualManifest.getFiles(), containsInAnyOrder(expectedManifest.getFiles().toArray()));
  }

}
