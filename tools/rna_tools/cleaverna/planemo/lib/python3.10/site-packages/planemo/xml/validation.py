"""Module describing abstractions for validating XML content."""

import abc
import subprocess
from collections import namedtuple

from galaxy.util import unicodify
from galaxy.util.commands import which

try:
    from lxml import etree
except ImportError:
    etree = None  # type: ignore

XMLLINT_COMMAND = "xmllint --noout --schema {0} {1} 2>&1"
INSTALL_VALIDATOR_MESSAGE = (
    "This feature requires an external dependency "
    "to function, pleaes install xmllint (e.g 'brew "
    "install libxml2' or 'apt-get install "
    "libxml2-utils'."
)


class XsdValidator(metaclass=abc.ABCMeta):
    """Class allowing validation of XML files against XSD schema."""

    @abc.abstractmethod
    def validate(self, schema_path, target_path):
        """Validate ``target_path`` against ``schema_path``.

        :return type: ValidationResult
        """

    @abc.abstractmethod
    def enabled(self):
        """Return True iff system has dependencies for this validator.

        :return type: bool
        """


ValidationResult = namedtuple("ValidationResult", ["passed", "output"])


class LxmlValidator(XsdValidator):
    """Validate XSD files using lxml library."""

    def validate(self, schema_path, target_path):
        try:
            xsd_doc = etree.parse(schema_path)
            xsd = etree.XMLSchema(xsd_doc)
            xml = etree.parse(target_path)
            passed = xsd.validate(xml)
            return ValidationResult(passed, xsd.error_log)
        except etree.XMLSyntaxError as e:
            return ValidationResult(False, unicodify(e))

    def enabled(self):
        return etree is not None


class XmllintValidator(XsdValidator):
    """Validate XSD files with the external tool xmllint."""

    def validate(self, schema_path, target_path):
        command = XMLLINT_COMMAND.format(schema_path, target_path)
        p = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
        stdout, _ = p.communicate()
        passed = p.returncode == 0
        return ValidationResult(passed, stdout)

    def enabled(self):
        return bool(which("xmllint"))


VALIDATORS = [LxmlValidator(), XmllintValidator()]


def get_validator(require=True):
    """Return a :class:`XsdValidator` object based on available dependencies."""
    for validator in VALIDATORS:
        if validator.enabled():
            return validator

    if require:
        raise Exception(INSTALL_VALIDATOR_MESSAGE)

    return None


__all__ = (
    "get_validator",
    "XsdValidator",
)
