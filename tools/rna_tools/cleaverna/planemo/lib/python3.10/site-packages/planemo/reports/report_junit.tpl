<?xml version="1.0" encoding="UTF-8"?>
<testsuites>
    <testsuite name="{{ raw_data.suitename }}"
              tests="{{ raw_data.results.total }}"
              errors="{{ raw_data.results.errors }}"
              failures="{{ raw_data.results.failures }}"
              skipped="{{ raw_data.results.skips }}">
        {% for testcase in raw_data.tests %}
        <testcase classname="{{ testcase.id }}" name="{{ testcase.data.test_index }}" time="{{ testcase.data.time_seconds }}">
            {% if 'job' in testcase.data %}
                {% if testcase.data.status != 'success' %}
                    <failure message="Tool exit code: {{ testcase.data.job.exit_code }}"><![CDATA[
                        {{ testcase.data | tojson(indent=True)| strip_control_characters }}
                    ]]></failure>
                {% endif %}
            {% else %}
                <failure message="{{  testcase.data.execution_problem }}"><![CDATA[
                    {{ testcase.data | tojson(indent=True)| strip_control_characters }}
                ]]></failure>
            {% endif %}
        </testcase>
        {% endfor %}
   </testsuite>
</testsuites>
