workflow Workflow {
    call makeSomeFiles {}
    call takesSomeFiles {input: files = makeSomeFiles.createdFiles}
    output {
        File outFile = takesSomeFiles.outFile
    }
}

task makeSomeFiles {
    command {
        for I in $(seq 1 25)
        do
          echo $I > some_filename_that_is_way_longer_than_it_has_to_be_because_that_may_or_may_not_be_relevant_to_triggering_the_behavior_this_test_is_supposed_to_be_testing_$I.txt
        done
    }

    output {
        Array[File] createdFiles = [
            "some_filename_that_is_way_longer_than_it_has_to_be_because_that_may_or_may_not_be_relevant_to_triggering_the_behavior_this_test_is_supposed_to_be_testing_1.txt",
            "some_filename_that_is_way_longer_than_it_has_to_be_because_that_may_or_may_not_be_relevant_to_triggering_the_behavior_this_test_is_supposed_to_be_testing_2.txt",
            "some_filename_that_is_way_longer_than_it_has_to_be_because_that_may_or_may_not_be_relevant_to_triggering_the_behavior_this_test_is_supposed_to_be_testing_3.txt",
            "some_filename_that_is_way_longer_than_it_has_to_be_because_that_may_or_may_not_be_relevant_to_triggering_the_behavior_this_test_is_supposed_to_be_testing_4.txt",
            "some_filename_that_is_way_longer_than_it_has_to_be_because_that_may_or_may_not_be_relevant_to_triggering_the_behavior_this_test_is_supposed_to_be_testing_5.txt",
            "some_filename_that_is_way_longer_than_it_has_to_be_because_that_may_or_may_not_be_relevant_to_triggering_the_behavior_this_test_is_supposed_to_be_testing_6.txt",
            "some_filename_that_is_way_longer_than_it_has_to_be_because_that_may_or_may_not_be_relevant_to_triggering_the_behavior_this_test_is_supposed_to_be_testing_7.txt",
            "some_filename_that_is_way_longer_than_it_has_to_be_because_that_may_or_may_not_be_relevant_to_triggering_the_behavior_this_test_is_supposed_to_be_testing_8.txt",
            "some_filename_that_is_way_longer_than_it_has_to_be_because_that_may_or_may_not_be_relevant_to_triggering_the_behavior_this_test_is_supposed_to_be_testing_9.txt",
            "some_filename_that_is_way_longer_than_it_has_to_be_because_that_may_or_may_not_be_relevant_to_triggering_the_behavior_this_test_is_supposed_to_be_testing_10.txt",
            "some_filename_that_is_way_longer_than_it_has_to_be_because_that_may_or_may_not_be_relevant_to_triggering_the_behavior_this_test_is_supposed_to_be_testing_11.txt",
            "some_filename_that_is_way_longer_than_it_has_to_be_because_that_may_or_may_not_be_relevant_to_triggering_the_behavior_this_test_is_supposed_to_be_testing_12.txt",
            "some_filename_that_is_way_longer_than_it_has_to_be_because_that_may_or_may_not_be_relevant_to_triggering_the_behavior_this_test_is_supposed_to_be_testing_13.txt",
            "some_filename_that_is_way_longer_than_it_has_to_be_because_that_may_or_may_not_be_relevant_to_triggering_the_behavior_this_test_is_supposed_to_be_testing_14.txt",
            "some_filename_that_is_way_longer_than_it_has_to_be_because_that_may_or_may_not_be_relevant_to_triggering_the_behavior_this_test_is_supposed_to_be_testing_15.txt",
            "some_filename_that_is_way_longer_than_it_has_to_be_because_that_may_or_may_not_be_relevant_to_triggering_the_behavior_this_test_is_supposed_to_be_testing_16.txt",
            "some_filename_that_is_way_longer_than_it_has_to_be_because_that_may_or_may_not_be_relevant_to_triggering_the_behavior_this_test_is_supposed_to_be_testing_17.txt",
            "some_filename_that_is_way_longer_than_it_has_to_be_because_that_may_or_may_not_be_relevant_to_triggering_the_behavior_this_test_is_supposed_to_be_testing_18.txt",
            "some_filename_that_is_way_longer_than_it_has_to_be_because_that_may_or_may_not_be_relevant_to_triggering_the_behavior_this_test_is_supposed_to_be_testing_19.txt",
            "some_filename_that_is_way_longer_than_it_has_to_be_because_that_may_or_may_not_be_relevant_to_triggering_the_behavior_this_test_is_supposed_to_be_testing_20.txt",
            "some_filename_that_is_way_longer_than_it_has_to_be_because_that_may_or_may_not_be_relevant_to_triggering_the_behavior_this_test_is_supposed_to_be_testing_21.txt",
            "some_filename_that_is_way_longer_than_it_has_to_be_because_that_may_or_may_not_be_relevant_to_triggering_the_behavior_this_test_is_supposed_to_be_testing_22.txt",
            "some_filename_that_is_way_longer_than_it_has_to_be_because_that_may_or_may_not_be_relevant_to_triggering_the_behavior_this_test_is_supposed_to_be_testing_23.txt",
            "some_filename_that_is_way_longer_than_it_has_to_be_because_that_may_or_may_not_be_relevant_to_triggering_the_behavior_this_test_is_supposed_to_be_testing_24.txt",
            "some_filename_that_is_way_longer_than_it_has_to_be_because_that_may_or_may_not_be_relevant_to_triggering_the_behavior_this_test_is_supposed_to_be_testing_25.txt"
        ]
    }
    runtime {
        docker: "ubuntu"
        memory: "2G"
        cpu: 1
  }
}

task takesSomeFiles {
    String outputFileName = "output.txt"
    Array[File] files

    command {
        echo 'Hello world' > ${outputFileName}
    }

    runtime {
        docker: "ubuntu"
        memory: "2G"
        cpu: 1
    }

    output {
        File outFile = outputFileName
    }
}

