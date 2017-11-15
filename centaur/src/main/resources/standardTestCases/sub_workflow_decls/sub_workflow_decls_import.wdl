task first_task {
  command {
    echo foo
  }
  runtime {
    docker: "ubuntu:latest"
  }
  output {
    String out = read_string(stdout())
  }
}

task second_task {
  String relay
  command {
    echo ${relay}
  }
  runtime {
    docker: "ubuntu:latest"
  }
  output {
    String out = read_string(stdout())
  }
}

workflow sub_decls {
  # A subworkflow input that needs to be supplied.
  String needs_to_be_supplied

  # A subworkflow input with a default which is overridden.
  String default_which_is_overridden = "don't want to see this"

  # A subworkflow declaration with a default that depends on another input which has no default.
  String depends_on_needs_to_be_supplied = needs_to_be_supplied

  # A subworkflow declaration with a default that depends on another input which is overridden.
  String depends_on_default_which_is_overridden = default_which_is_overridden

  # An overridden subworkflow declaration with a default that depends on another input which has no default.
  String overridden_depends_on_needs_to_be_supplied = needs_to_be_supplied

  # An overridden subworkflow declaration with a default that depends on another input which was overridden.
  String overridden_depends_on_default_which_is_overridden = default_which_is_overridden

  # An optional subworkflow input no default
  # Uninitialized optionals in subworkflows currently crash and burn (the SubWorkflowExecutionActor dies
  # without diagnostics).
  # String? optional_without_default_not_supplied
  # String? optional_without_default_supplied

  # An optional subworkflow input that's got a default
  String? optional_with_default = "this is a default"

  # An optional subworkflow input that's got a default but is overridden
  String? optional_with_default_but_overridden = "don't want to see this optional"

  # Passthrough output from the outer workflow.
  String passthrough_required_input = "Don't want to see this value in the output!!!"
  String passthrough_initialized_input = "Don't want to see this value in the output!!!"
  String passthrough_overridden_input = "Don't want to see this value in the output!!!"

  call first_task
  # This subworkflow declaration is initialized and not overridden, so it shouldn't be
  # evaluated at subworkflow call time. If it is evaluated at subworkflow call time that
  # will fail since the task call on which this depends hasn't run yet.
  # This is the most recent fix to the 29 hotfix branch.
  String depends_on_task_output = first_task.out
  call second_task { input: relay = depends_on_task_output }
  output {
    String depends_on_task_output_out = second_task.out
    String needs_to_be_supplied_out = needs_to_be_supplied
    String default_which_is_overridden_out = default_which_is_overridden
    String depends_on_needs_to_be_supplied_out = depends_on_needs_to_be_supplied
    String depends_on_default_which_is_overridden_out = depends_on_default_which_is_overridden
    String overridden_depends_on_needs_to_be_supplied_out = overridden_depends_on_needs_to_be_supplied
    String overridden_depends_on_default_which_is_overridden_out = overridden_depends_on_default_which_is_overridden
    # String? optional_without_default_not_supplied_out = optional_without_default_not_supplied
    # String? optional_without_default_supplied_out = optional_without_default_supplied
    String? optional_with_default_out = optional_with_default
    String? optional_with_default_but_overridden_out = optional_with_default_but_overridden
    String passthrough_required_input_out = passthrough_required_input
    String passthrough_initialized_input_out = passthrough_initialized_input
    String passthrough_overridden_input_out = passthrough_overridden_input
  }
}
