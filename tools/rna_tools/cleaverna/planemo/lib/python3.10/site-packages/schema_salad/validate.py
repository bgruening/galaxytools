import logging
import pprint
from collections.abc import Mapping, MutableMapping, MutableSequence
from typing import Any, Final, Optional
from urllib.parse import urlsplit

from . import avro
from .avro.schema import Schema
from .exceptions import (
    ClassValidationException,
    SchemaSaladException,
    ValidationException,
)
from .sourceline import SourceLine

_logger: Final = logging.getLogger("salad")


def validate(
    expected_schema: Schema,
    datum: Any,
    identifiers: Optional[list[str]] = None,
    strict: bool = False,
    foreign_properties: Optional[set[str]] = None,
    vocab: Optional[Mapping[str, str]] = None,
) -> bool:
    if not identifiers:
        identifiers = []
    if not foreign_properties:
        foreign_properties = set()
    return validate_ex(
        expected_schema,
        datum,
        identifiers,
        strict=strict,
        foreign_properties=foreign_properties,
        raise_ex=False,
        vocab=vocab,
    )


INT_MIN_VALUE: Final = -(1 << 31)
INT_MAX_VALUE: Final = (1 << 31) - 1
LONG_MIN_VALUE: Final = -(1 << 63)
LONG_MAX_VALUE: Final = (1 << 63) - 1


def avro_shortname(name: str) -> str:
    """Produce an avro friendly short name."""
    return name.split(".")[-1]


saladp: Final = "https://w3id.org/cwl/salad#"
primitives: Final = {
    "http://www.w3.org/2001/XMLSchema#string": "string",
    "http://www.w3.org/2001/XMLSchema#boolean": "boolean",
    "http://www.w3.org/2001/XMLSchema#int": "int",
    "http://www.w3.org/2001/XMLSchema#long": "long",
    "http://www.w3.org/2001/XMLSchema#float": "float",
    "http://www.w3.org/2001/XMLSchema#double": "double",
    saladp + "null": "null",
    saladp + "enum": "enum",
    saladp + "array": "array",
    saladp + "record": "record",
    saladp + "map": "map",
    saladp + "union": "union",
}


def avro_type_name(url: str) -> str:
    """
    Turn a URL into an Avro-safe name.

    If the URL has no fragment, return this plain URL.

    Extract either the last part of the URL fragment past the slash, otherwise
    the whole fragment.
    """

    if url in primitives:
        return primitives[url]

    u: Final = urlsplit(url)
    joined: Final = filter(
        lambda x: x,
        list(reversed(u.netloc.split("."))) + u.path.split("/") + u.fragment.split("/"),
    )
    return ".".join(joined)


def friendly(v: Any) -> Any:
    """Format an Avro schema into a pretty-printed representation."""
    if isinstance(v, avro.schema.NamedSchema):
        return avro_shortname(v.name)
    if isinstance(v, avro.schema.ArraySchema):
        return f"array of <{friendly(v.items)}>"
    if isinstance(v, (avro.schema.MapSchema, avro.schema.NamedMapSchema)):
        return f"map of <{friendly(v.values)}>"
    if isinstance(v, avro.schema.PrimitiveSchema):
        return v.type
    if isinstance(v, (avro.schema.UnionSchema, avro.schema.NamedUnionSchema)):
        return " or ".join([friendly(s) for s in v.schemas])
    return avro_shortname(v)


def vpformat(datum: Any) -> str:
    """Truncate a pretty-printed representation of a Python object to 160 characters."""
    a = pprint.pformat(datum)
    if len(a) > 160:
        a = a[0:160] + "[...]"
    return a


def validate_ex(
    expected_schema: Schema,
    datum: Any,
    identifiers: Optional[list[str]] = None,
    strict: bool = False,
    foreign_properties: Optional[set[str]] = None,
    raise_ex: bool = True,
    strict_foreign_properties: bool = False,
    logger: logging.Logger = _logger,
    skip_foreign_properties: bool = False,
    vocab: Optional[Mapping[str, str]] = None,
) -> bool:
    """Determine if a python datum is an instance of a schema."""
    debug: Final = _logger.isEnabledFor(logging.DEBUG)
    if not identifiers:
        identifiers = []

    if not foreign_properties:
        foreign_properties = set()

    if vocab is None:
        raise Exception("vocab must be provided")

    schema_type: Final = expected_schema.type

    if schema_type == "null":
        if datum is None:
            return True
        if raise_ex:
            raise ValidationException("the value is not null")
        return False
    if schema_type == "boolean":
        if isinstance(datum, bool):
            return True
        if raise_ex:
            raise ValidationException("the value is not boolean")
        return False
    if schema_type == "string":
        if isinstance(datum, str):
            return True
        if isinstance(datum, bytes):
            return True
        if raise_ex:
            raise ValidationException("the value is not string")
        return False
    if schema_type == "int":
        if isinstance(datum, int) and INT_MIN_VALUE <= datum <= INT_MAX_VALUE:
            return True
        if raise_ex:
            raise ValidationException(f"{vpformat(datum)!r} is not int")
        return False
    if schema_type == "long":
        if (isinstance(datum, int)) and LONG_MIN_VALUE <= datum <= LONG_MAX_VALUE:
            return True
        if raise_ex:
            raise ValidationException(f"the value {vpformat(datum)!r} is not long")
        return False
    if schema_type in ["float", "double"]:
        if isinstance(datum, (int, float)):
            return True
        if raise_ex:
            raise ValidationException(f"the value {vpformat(datum)!r} is not float or double")
        return False
    if isinstance(expected_schema, avro.schema.EnumSchema):
        if expected_schema.name in ("org.w3id.cwl.salad.Any", "Any"):
            if datum is not None:
                return True
            if raise_ex:
                raise ValidationException("'Any' type must be non-null")
            return False
        if not isinstance(datum, str):
            if raise_ex:
                raise ValidationException(
                    f"value is a {type(datum).__name__} but expected a string"
                )
            return False
        if expected_schema.name == "org.w3id.cwl.cwl.Expression":
            if "$(" in datum or "${" in datum:
                return True
            if raise_ex:
                raise ValidationException(
                    f"value {datum!r} does not contain an expression in the form $() or ${{}}"
                )
            return False
        if datum in expected_schema.symbols:
            return True
        if raise_ex:
            raise ValidationException(
                "the value {} is not a valid {}, expected {}{}".format(
                    vpformat(datum),
                    friendly(expected_schema.name),
                    "one of " if len(expected_schema.symbols) > 1 else "",
                    "'" + "', '".join(expected_schema.symbols) + "'",
                )
            )
        return False
    if isinstance(expected_schema, avro.schema.ArraySchema):
        if isinstance(datum, MutableSequence):
            for i, d in enumerate(datum):
                try:
                    sl = SourceLine(datum, i, ValidationException)
                    if not validate_ex(
                        expected_schema.items,
                        d,
                        identifiers,
                        strict=strict,
                        foreign_properties=foreign_properties,
                        raise_ex=raise_ex,
                        strict_foreign_properties=strict_foreign_properties,
                        logger=logger,
                        skip_foreign_properties=skip_foreign_properties,
                        vocab=vocab,
                    ):
                        return False
                except ValidationException as v:
                    if raise_ex:
                        source = v if debug else None
                        raise ValidationException("item is invalid because", sl, [v]) from source
                    return False
            return True
        if raise_ex:
            raise ValidationException(
                f"the value {vpformat(datum)} is not a list, "
                f"expected list of {friendly(expected_schema.items)}"
            )
        return False
    if isinstance(expected_schema, (avro.schema.UnionSchema, avro.schema.NamedUnionSchema)):
        for s in expected_schema.schemas:
            if validate_ex(
                s,
                datum,
                identifiers,
                strict=strict,
                raise_ex=False,
                strict_foreign_properties=strict_foreign_properties,
                logger=logger,
                skip_foreign_properties=skip_foreign_properties,
                vocab=vocab,
            ):
                return True

        if not raise_ex:
            return False

        errors1: Final[list[SchemaSaladException]] = []
        checked: Final = []
        for s in expected_schema.schemas:
            if isinstance(datum, MutableSequence) and not isinstance(s, avro.schema.ArraySchema):
                continue
            if isinstance(datum, MutableMapping) and not isinstance(
                s, (avro.schema.RecordSchema, avro.schema.MapSchema, avro.schema.NamedMapSchema)
            ):
                continue
            if isinstance(datum, (bool, int, float, str)) and isinstance(
                s,
                (
                    avro.schema.ArraySchema,
                    avro.schema.RecordSchema,
                    avro.schema.MapSchema,
                    avro.schema.NamedMapSchema,
                ),
            ):
                continue
            if datum is not None and s.type == "null":
                continue

            checked.append(s)
            try:
                validate_ex(
                    s,
                    datum,
                    identifiers,
                    strict=strict,
                    foreign_properties=foreign_properties,
                    raise_ex=True,
                    strict_foreign_properties=strict_foreign_properties,
                    logger=logger,
                    skip_foreign_properties=skip_foreign_properties,
                    vocab=vocab,
                )
            except ClassValidationException:
                raise
            except ValidationException as e:
                errors1.append(e)
        if bool(errors1):
            raise ValidationException(
                "",
                None,
                [
                    ValidationException(f"tried {friendly(check)} but", None, [err])
                    for (check, err) in zip(checked, errors1)
                ],
                "-",
            )
        raise ValidationException(
            f"value is a {type(datum).__name__}, expected {friendly(expected_schema)}"
        )

    if isinstance(expected_schema, avro.schema.RecordSchema):
        if not isinstance(datum, MutableMapping):
            if raise_ex:
                raise ValidationException(
                    f"is not a dict. Expected a {friendly(expected_schema.name)} object."
                )
            return False

        classmatch = None
        for f in expected_schema.fields:
            if f.name in ("class",):
                d = datum.get(f.name)
                if not d:
                    if raise_ex:
                        raise ValidationException(f"Missing {f.name!r} field")
                    return False
                avroname = None
                if d in vocab:
                    avroname = avro_type_name(vocab[d])
                if expected_schema.name not in (d, avroname):
                    if raise_ex:
                        raise ValidationException(
                            f"Expected class {expected_schema.name!r} but this is {d!r}"
                        )
                    return False
                classmatch = d
                break

        errors2: Final = []
        for f in expected_schema.fields:
            if f.name in ("class",):
                continue

            if f.name in datum:
                fieldval = datum[f.name]
            else:
                try:
                    fieldval = f.default
                except KeyError:
                    fieldval = None

            try:
                sl = SourceLine(datum, f.name, str)
                if not validate_ex(
                    f.type,
                    fieldval,
                    identifiers,
                    strict=strict,
                    foreign_properties=foreign_properties,
                    raise_ex=raise_ex,
                    strict_foreign_properties=strict_foreign_properties,
                    logger=logger,
                    skip_foreign_properties=skip_foreign_properties,
                    vocab=vocab,
                ):
                    return False
            except ValidationException as v:
                if f.name not in datum:
                    errors2.append(ValidationException(f"missing required field {f.name!r}"))
                else:
                    errors2.append(
                        ValidationException(
                            f"the {f.name!r} field is not valid because",
                            sl,
                            [v],
                        )
                    )

        for d in datum:
            found = False
            for f in expected_schema.fields:
                if d == f.name:
                    found = True
            if not found:
                sl = SourceLine(datum, d, str)
                if d is None:
                    err = ValidationException("mapping with implicit null key", sl)
                    if strict:
                        errors2.append(err)
                    else:
                        logger.warning(err.as_warning())
                    continue
                if d not in identifiers and d not in foreign_properties and d[0] not in ("@", "$"):
                    if (
                        (d not in identifiers and strict)
                        and (
                            d not in foreign_properties
                            and strict_foreign_properties
                            and not skip_foreign_properties
                        )
                        and not raise_ex
                    ):
                        return False
                    split = urlsplit(d)
                    if split.scheme:
                        if not skip_foreign_properties:
                            err = ValidationException(
                                "unrecognized extension field {!r}{}.{}".format(
                                    d,
                                    (
                                        " and strict_foreign_properties checking is enabled"
                                        if strict_foreign_properties
                                        else ""
                                    ),
                                    (
                                        "\nForeign properties from $schemas:\n  {}".format(
                                            "\n  ".join(sorted(foreign_properties))
                                        )
                                        if len(foreign_properties) > 0
                                        else ""
                                    ),
                                ),
                                sl,
                            )
                            if strict_foreign_properties:
                                errors2.append(err)
                            elif len(foreign_properties) > 0:
                                logger.warning(err.as_warning())
                    else:
                        err = ValidationException(
                            "invalid field {!r}, expected one of: {}".format(
                                d,
                                ", ".join(f"{fn.name!r}" for fn in expected_schema.fields),
                            ),
                            sl,
                        )
                        if strict:
                            errors2.append(err)
                        else:
                            logger.warning(err.as_warning())

        if bool(errors2):
            if raise_ex:
                if classmatch:
                    raise ClassValidationException("", None, errors2, "*")
                raise ValidationException("", None, errors2, "*")
            return False
        return True

    if isinstance(expected_schema, (avro.schema.MapSchema, avro.schema.NamedMapSchema)):
        if isinstance(datum, MutableMapping):
            for key, val in datum.items():
                if not isinstance(key, str):
                    pass
                try:
                    sl = SourceLine(datum, key, ValidationException)
                    if not validate_ex(
                        expected_schema.values,
                        val,
                        identifiers,
                        strict=strict,
                        foreign_properties=foreign_properties,
                        raise_ex=raise_ex,
                        strict_foreign_properties=strict_foreign_properties,
                        logger=logger,
                        skip_foreign_properties=skip_foreign_properties,
                        vocab=vocab,
                    ):
                        return False
                except ValidationException as v:
                    if raise_ex:
                        source = v if debug else None
                        raise ValidationException("item is invalid because", sl, [v]) from source
                    return False
            return True
        if raise_ex:
            raise ValidationException(
                f"the value {vpformat(datum)} is not an object, "
                f"expected object of {friendly(expected_schema.values)}"
            )
        return False

    if raise_ex:
        raise ValidationException(f"Unrecognized schema_type {schema_type}")
    return False
