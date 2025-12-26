{% from 'macros.tmpl' import render_invocation_details, render_invocation_messages, render_job_parameters, render_steps %}
{% if title %}
# {{ execution_type }} {{ title }}

{% endif %}
## {{ execution_type }} Summary
{% set state = namespace(found=false) %}
{% set state.success = raw_data.results.total - raw_data.results.errors - raw_data.results.failures - raw_data.results.skips | default(0) %}
{% set state.error = raw_data.results.errors | default(0) %}
{% set state.failure = raw_data.results.failures | default(0) %}
{% set state.skipped = raw_data.results.skipped | default(0) %}

{% if raw_data.results.total %}
<div class="progress">
  <div class="progress-bar progress-bar-success" style="width: {{ (state.success / raw_data.results.total) * 100 }}%" aria-valuenow="{{ state.success }}" aria-valuemin="0" aria-valuemax="{{ raw_data.results.total }}" data-toggle="tooltip" title="{{state.success}} Passed">
  </div>
  <div class="progress-bar progress-bar-warning" style="width: {{ (state.skipped / raw_data.results.total) * 100 }}%" aria-valuenow="{{ state.skipped }}" aria-valuemin="0" aria-valuemax="{{ raw_data.results.total }}" data-toggle="tooltip" title="{{state.skipped}} Skipped">
  </div>
  <div class="progress-bar progress-bar-danger" style="width: {{ ((state.error + state.failure) / raw_data.results.total) * 100 }}%" aria-valuenow="{{ state.error + state.failure }}" aria-valuemin="0" aria-valuemax="{{ raw_data.results.total }}" title="{{state.error + state.failure}} Failed or Errored">
  </div>
</div>
{% endif %}

| {{ execution_type }} State | Count |
| ---------- | ----- |
| Total      | {{ raw_data.results.total | default(0)  }} |
| Passed     | {{ state.success }} |
| Error      | {{ state.error }} |
| Failure    | {{ state.failure }} |
| Skipped    | {{ state.skipped }} |

{% for status, desc in {'error': 'Errored', 'failure': 'Failed', 'success': 'Passed'}.items() if state[status]%}
{% set expanded = "open" if status in ("error", "failure") else "" %}
<details {{ expanded }}><summary>{{ desc }} {{ execution_type }}s</summary>
{%   for test in raw_data.tests %}
{%     if test.data.status == status %}
{%       if test.data.status == 'success' %}
* <summary class="light-green">&#9989; {{ test.id|replace("#","# ") }}</summary>
{%       else %}
* <summary class="light-red">&#10060; {{ test.id|replace("#","# ") }}</summary>
{%       endif %}
{%     endif %}
{%   endfor %}

</details>
{% endfor %}
