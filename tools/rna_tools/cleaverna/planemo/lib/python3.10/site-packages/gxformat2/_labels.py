"""Utilities for handling unlabelled objects when translating workflow formats."""


class Labels:
    """Track labels assigned and generate anonymous ones."""

    def __init__(self):
        """Initialize labels that have been encountered or generated."""
        self.seen_labels = set()
        self.anonymous_labels = 0

    def ensure_new_output_label(self, label: str):
        """Ensure supplied label has value or generate an anonymous one."""
        if label is None:
            self.anonymous_labels += 1
            label = f"_anonymous_output_{self.anonymous_labels}"
        assert label not in self.seen_labels
        self.seen_labels.add(label)
        return label

    @staticmethod
    def is_anonymous_output_label(label: str):
        """Predicate determining if supplied label was generated for anonymous output."""
        # label likely can't be null according to the schema definition - but we've got a test
        # in Galaxy that doesn't define a label in order to great a .ga file without output
        # labels (which is completely normal).
        return not label or label.startswith("_anonymous_output_")
