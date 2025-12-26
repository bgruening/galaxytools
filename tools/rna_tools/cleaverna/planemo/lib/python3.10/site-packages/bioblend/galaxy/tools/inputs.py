from collections.abc import Iterator
from typing import (
    Any,
    Optional,
    Union,
)


class InputsBuilder:
    """ """

    def __init__(self) -> None:
        self._input_dict: dict[str, Any] = {}

    def set(self, name: str, input: Any) -> "InputsBuilder":
        self._input_dict[name] = input
        return self

    def set_param(self, name: str, value: Any) -> "InputsBuilder":
        return self.set(name, param(value=value))

    def set_dataset_param(self, name: str, value: str, src: str = "hda") -> "InputsBuilder":
        return self.set(name, dataset(value, src=src))

    def to_dict(self) -> dict[str, Any]:
        values = {}
        for key, value in self.flat_iter():
            if hasattr(value, "value"):
                value = value.value
            values[key] = value
        return values

    def flat_iter(self, prefix: Optional[str] = None) -> Iterator[tuple[str, Any]]:
        for key, value in self._input_dict.items():
            effective_key = key if prefix is None else f"{prefix}|{key}"
            if hasattr(value, "flat_iter"):
                yield from value.flat_iter(effective_key)
            else:
                yield effective_key, value


class RepeatBuilder:
    def __init__(self) -> None:
        self._instances: list[InputsBuilder] = []

    def instance(self, inputs: InputsBuilder) -> "RepeatBuilder":
        self._instances.append(inputs)
        return self

    def flat_iter(self, prefix: str) -> Iterator[tuple[str, Any]]:
        for index, instance in enumerate(self._instances):
            index_prefix = f"{prefix}_{index}"
            yield from instance.flat_iter(index_prefix)


class Param:
    def __init__(self, value: Any) -> None:
        self.value = value


class DatasetParam(Param):
    def __init__(self, value: Union[dict[str, str], str], src: str = "hda") -> None:
        if not isinstance(value, dict):
            value = {"src": src, "id": value}
        super().__init__(value)


inputs = InputsBuilder
repeat = RepeatBuilder
conditional = InputsBuilder
param = Param
dataset = DatasetParam

__all__ = ("inputs", "repeat", "conditional", "param")
