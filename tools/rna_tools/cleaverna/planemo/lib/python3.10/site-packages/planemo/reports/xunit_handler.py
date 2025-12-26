from .build_report import template_data


def handle_report_xunit_kwd(kwds, collected_data):
    if kwds.get("report_xunit", False):
        with open(kwds["report_xunit"], "w") as handle:
            handle.write(template_data(collected_data, template_name="xunit.tpl"))
