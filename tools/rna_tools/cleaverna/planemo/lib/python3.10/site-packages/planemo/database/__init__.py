"""Package is responsible for managing planemo profile databases.

This package makes it very easy to create and destroy databases, therefore it
should not be used for production data - and should not even be connnected
to a production database server.
"""

from .factory import create_database_source

__all__ = ("create_database_source",)
