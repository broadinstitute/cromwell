version development

workflow test_custom_cacheworthy_attributes {
    call custom_cacheworthy_attributes_task_origin
    call custom_cacheworthy_attributes_task_should_cache after custom_cacheworthy_attributes_task_origin
    call custom_cacheworthy_attributes_task_should_not_cache after custom_cacheworthy_attributes_task_origin
}

task custom_cacheworthy_attributes_task_origin {
    input { }
    runtime {
        cacheworthy_attribute: 15
        uncacheworthy_attribute_1: 16
        uncacheworthy_attribute_2: 17
        # Config for this backend is defined in src/ci/resources/local_application.conf
        backend: "LocalCacheableRuntimeAttribute"
    }
    command {
      echo foo > foo
    }
    output {
        String foo = read_string("foo")
    }
}

# This task should cache to the 'origin' even though the two "uncacheable" attributes are different
task custom_cacheworthy_attributes_task_should_cache {
    input { }
    runtime {
        cacheworthy_attribute: 15
        uncacheworthy_attribute_1: 160
        uncacheworthy_attribute_2: 170
        # Config for this backend is defined in src/ci/resources/local_application.conf
        backend: "LocalCacheableRuntimeAttribute"
    }
    command {
      echo foo > foo
    }
    output {
        String foo = read_string("foo")
    }
}

# This task should not cache to the 'origin' because the "cacheable" attribute is different
task custom_cacheworthy_attributes_task_should_not_cache {
    input { }
    runtime {
        cacheworthy_attribute: 150
        uncacheworthy_attribute_1: 16
        uncacheworthy_attribute_2: 17
        # Config for this backend is defined in src/ci/resources/local_application.conf
        backend: "LocalCacheableRuntimeAttribute"
    }
    command {
      echo foo > foo
    }
    output {
        String foo = read_string("foo")
    }
}
