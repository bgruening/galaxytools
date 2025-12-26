"""Describe the interface classes of the planemo.database package."""

import abc


class DatabaseSource(metaclass=abc.ABCMeta):
    """Interface describing a source of profile databases."""

    @abc.abstractmethod
    def create_database(self, identifier):
        """Create a database with specified short identifier.

        Throw an exception if it already exists.
        """

    @abc.abstractmethod
    def delete_database(self, identifier):
        """Delete a database with specified short identifier.

        Throw an exception if it already exists.
        """

    @abc.abstractmethod
    def list_databases(self):
        """Return identifiers associated with database source."""

    @abc.abstractmethod
    def sqlalchemy_url(self, identifier):
        """Return a URL string for use by sqlalchemy."""


__all__ = ("DatabaseSource",)
