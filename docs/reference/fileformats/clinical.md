::: genie_registry.clinical
    options:
      members:
        - _check_year
        - _check_int_dead_consistency
        - _check_int_year_consistency
        - _check_year_death_validity
        - _check_year_death_validity_message
        - _check_int_dod_validity
        - _check_int_dod_validity_message
        - remap_clinical_values
      members_order: source
      show_if_no_docstring: true

::: genie_registry.clinical.Clinical
    options:
      members: true
      members_order: source
      show_if_no_docstring: true
