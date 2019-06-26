version 1.0

workflow wf_level_file_size {
    File input1 = "dos://wb-mock-drs-dev.storage.googleapis.com/4a3908ad-1f0b-4e2a-8a92-611f2123e8b0"
    File input2 = "dos://wb-mock-drs-dev.storage.googleapis.com/0c8e7bc6-fd76-459d-947b-808b0605beb3"

    output {
        Float fileSize1 = size(input1)
        Float fileSize2 = size(input2)
    }
}
