{% set bad = raw_data.summary.num_errors + raw_data.summary.num_failures + raw_data.summary.num_skips -%}
{% if bad == 0 -%}
All {{ raw_data.summary.num_tests }} test(s) executed passed.
{%- else -%}
There were problems with {{ bad }} out of {{raw_data.summary.num_tests}} test(s) executed.
{%- endif %}

{% for test in raw_data.tests -%}
{{ test.id | replace('functional.test_toolbox.TestForTool_', '') }}: {{ test.data.status }}
{% if test.data.output_problems -%}
{% for problem in test.data.output_problems -%}
{{problem | indent(2, True) }}
{% endfor %}
{% endif %}
{%- endfor %}
