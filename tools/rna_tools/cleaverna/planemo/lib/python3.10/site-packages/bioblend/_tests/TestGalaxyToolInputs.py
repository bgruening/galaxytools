from bioblend.galaxy.tools.inputs import (
    conditional,
    dataset,
    inputs,
    repeat,
)


def test_conditional():
    # Build up example inputs for random_lines1
    as_dict = (
        inputs()
        .set("num_lines", 5)
        .set("input", dataset("encoded1"))
        .set("seed_source", conditional().set("seed_source_selector", "set_seed").set("seed", "asdf"))
        .to_dict()
    )
    assert as_dict["num_lines"] == 5
    assert as_dict["input"]["src"] == "hda"
    assert as_dict["input"]["id"] == "encoded1"
    assert as_dict["seed_source|seed_source_selector"] == "set_seed"
    assert as_dict["seed_source|seed"] == "asdf"


def test_repeat():
    # Build up inputs for cat1
    as_dict = (
        inputs()
        .set("input1", dataset("encoded1"))
        .set(
            "queries",
            repeat()
            .instance(inputs().set_dataset_param("input2", "encoded2"))
            .instance(inputs().set_dataset_param("input2", "encoded3")),
        )
        .to_dict()
    )
    assert as_dict["input1"]["src"] == "hda"
    assert as_dict["input1"]["id"] == "encoded1"
    assert as_dict["queries_0|input2"]["src"] == "hda"
    assert as_dict["queries_0|input2"]["id"] == "encoded2"
    assert as_dict["queries_1|input2"]["src"] == "hda"
    assert as_dict["queries_1|input2"]["id"] == "encoded3"
