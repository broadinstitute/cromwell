import "sub_workflow_decls_import.wdl" as decls

workflow sub_workflow_decls {
  String workflow_required_input
  String workflow_initialized_input = "This is the correct initialized value"
  # This is somewhat stretching the definition of "overridden" as this optional is not initialized, but in 29
  # initialized optional workflow declarations cannot be overridden.
  String? workflow_overridden_input_optional
  String? workflow_default_input_optional = "Should not see this in the output"

  String workflow_overridden_input = select_first([workflow_overridden_input_optional, workflow_default_input_optional])

  call decls.sub_decls as sudecls { input:
    needs_to_be_supplied = "initialized"
    , default_which_is_overridden = "overridden"
    , overridden_depends_on_needs_to_be_supplied = "you've been overridden!"
    , overridden_depends_on_default_which_is_overridden = "your override has been overridden"
    #, optional_without_default_supplied = "consider yourself supplied"
    , optional_with_default_but_overridden = "supplied but overridden optional"
    , passthrough_required_input = workflow_required_input
    , passthrough_initialized_input = workflow_initialized_input
    , passthrough_overridden_input = workflow_overridden_input
  }

  output {
    String  depends_on_task_output = sudecls.depends_on_task_output_out
    String  needs_to_be_supplied = sudecls.needs_to_be_supplied_out
    String  default_which_is_overridden = sudecls.default_which_is_overridden_out
    String  depends_on_needs_to_be_supplied = sudecls.depends_on_needs_to_be_supplied_out
    String  depends_on_default_which_is_overridden = sudecls.depends_on_default_which_is_overridden_out
    String  overridden_depends_on_needs_to_be_supplied = sudecls.overridden_depends_on_needs_to_be_supplied_out
    String  overridden_depends_on_default_which_is_overridden = sudecls.overridden_depends_on_default_which_is_overridden_out
    # String? optional_without_default_not_supplied = sudecls.optional_without_default_not_supplied_out
    # String? optional_without_default_supplied = sudecls.optional_without_default_supplied_out
    String? optional_with_default = sudecls.optional_with_default_out
    String? optional_with_default_but_overridden = sudecls.optional_with_default_but_overridden_out
    String passthrough_required_input = sudecls.passthrough_required_input_out
    String passthrough_initialized_input = sudecls.passthrough_initialized_input_out
    String passthrough_overridden_input = sudecls.passthrough_overridden_input_out
  }
}
