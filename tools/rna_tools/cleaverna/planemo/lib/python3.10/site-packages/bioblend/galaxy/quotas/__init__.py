"""
Contains possible interactions with the Galaxy Quota
"""

from typing import (
    Any,
    Literal,
    Optional,
    TYPE_CHECKING,
)

from bioblend.galaxy.client import Client

if TYPE_CHECKING:
    from bioblend.galaxy import GalaxyInstance

QuotaOperations = Literal["+", "-", "="]
DefaultQuotaValues = Literal["no", "registered", "unregistered"]


class QuotaClient(Client):
    module = "quotas"

    def __init__(self, galaxy_instance: "GalaxyInstance") -> None:
        super().__init__(galaxy_instance)

    def get_quotas(self, deleted: bool = False) -> list[dict[str, Any]]:
        """
        Get a list of quotas

        :type deleted: bool
        :param deleted: Only return quota(s) that have been deleted

        :rtype: list
        :return: A list of dicts with details on individual quotas.
          For example::

            [{'id': '0604c8a56abe9a50',
              'model_class': 'Quota',
              'name': 'test ',
              'url': '/api/quotas/0604c8a56abe9a50'},
             {'id': '1ee267091d0190af',
              'model_class': 'Quota',
              'name': 'workshop',
              'url': '/api/quotas/1ee267091d0190af'}]
        """
        return self._get(deleted=deleted)

    def show_quota(self, quota_id: str, deleted: bool = False) -> dict[str, Any]:
        """
        Display information on a quota

        :type quota_id: str
        :param quota_id: Encoded quota ID

        :type deleted: bool
        :param deleted: Search for quota in list of ones already marked as deleted

        :rtype: dict
        :return: A description of quota.
          For example::

            {'bytes': 107374182400,
             'default': [],
             'description': 'just testing',
             'display_amount': '100.0 GB',
             'groups': [],
             'id': '0604c8a56abe9a50',
             'model_class': 'Quota',
             'name': 'test ',
             'operation': '=',
             'users': []}
        """
        return self._get(id=quota_id, deleted=deleted)

    def create_quota(
        self,
        name: str,
        description: str,
        amount: str,
        operation: QuotaOperations,
        default: Optional[DefaultQuotaValues] = "no",
        in_users: Optional[list[str]] = None,
        in_groups: Optional[list[str]] = None,
        quota_source_label: Optional[str] = None,
    ) -> dict[str, Any]:
        """
        Create a new quota

        :type name: str
        :param name: Name for the new quota. This must be unique within a Galaxy instance.

        :type description: str
        :param description: Quota description

        :type amount: str
        :param amount: Quota size (E.g. ``10000MB``, ``99 gb``, ``0.2T``, ``unlimited``)

        :type operation: str
        :param operation: One of (``+``, ``-``, ``=``)

        :type default: str
        :param default: Whether or not this is a default quota. Valid values
                        are "no", "unregistered", "registered" and None. None is
                        equivalent to "no".

        :type in_users: list of str
        :param in_users: A list of user IDs or user emails.

        :type in_groups: list of str
        :param in_groups: A list of group IDs or names.

        :type quota_source_label: str
        :param quota_source_label: If set, quota source label to apply this
          quota operation to. Otherwise, the default quota is used.

        :rtype: dict
        :return: A description of quota.
          For example::

            {'url': '/galaxy/api/quotas/386f14984287a0f7',
             'model_class': 'Quota',
             'message': "Quota 'Testing' has been created with 1 associated users and 0 associated groups.",
             'id': '386f14984287a0f7',
             'name': 'Testing'}
        """
        payload: dict[str, Any] = {
            "name": name,
            "description": description,
            "amount": amount,
            "operation": operation,
            "default": default,
            "quota_source_label": quota_source_label,
        }
        if in_users:
            payload["in_users"] = in_users

        if in_groups:
            payload["in_groups"] = in_groups

        return self._post(payload)

    def update_quota(
        self,
        quota_id: str,
        name: Optional[str] = None,
        description: Optional[str] = None,
        amount: Optional[str] = None,
        operation: Optional[QuotaOperations] = "=",
        default: Optional[DefaultQuotaValues] = None,
        in_users: Optional[list[str]] = None,
        in_groups: Optional[list[str]] = None,
    ) -> str:
        """
        Update an existing quota

        :type quota_id: str
        :param quota_id: Encoded quota ID

        :type name: str
        :param name: Name for the new quota. This must be unique within a Galaxy instance.

        :type description: str
        :param description: Quota description. If you supply this parameter,
                            but not the name, an error will be thrown.

        :type amount: str
        :param amount: Quota size (E.g. ``10000MB``, ``99 gb``, ``0.2T``, ``unlimited``)

        :type operation: str
        :param operation: One of (``+``, ``-``, ``=``). If you wish to change this
                          value, you must also provide the ``amount``,
                          otherwise it will not take effect.

        :type default: str
        :param default: Whether or not this is a default quota. Valid values
                        are "no", "unregistered", "registered" and None.
                        Calling this method with ``default="no"`` on a
                        non-default quota will throw an error. Passing None is
                        equivalent to not changing the current status.

        :type in_users: list of str
        :param in_users: A list of user IDs or user emails.

        :type in_groups: list of str
        :param in_groups: A list of group IDs or names.

        :rtype: str
        :return: A semicolon separated list of changes to the quota.
          For example::

            "Quota 'Testing-A' has been renamed to 'Testing-B'; Quota 'Testing-e' is now '-100.0 GB'; Quota 'Testing-B' is now the default for unregistered users"
        """
        payload: dict[str, Any] = {"default": default}
        if name:
            payload["name"] = name

        if description:
            payload["description"] = description

        if amount:
            payload["amount"] = amount

        if operation:
            payload["operation"] = operation

        if in_users:
            payload["in_users"] = in_users

        if in_groups:
            payload["in_groups"] = in_groups

        return self._put(id=quota_id, payload=payload)

    def delete_quota(self, quota_id: str) -> str:
        """
        Delete a quota

        Before a quota can be deleted, the quota must not be a default quota.

        :type quota_id: str
        :param quota_id: Encoded quota ID.

        :rtype: str
        :return: A description of the changes, mentioning the deleted quota.
          For example::

            "Deleted 1 quotas: Testing-B"
        """
        return self._delete(id=quota_id)

    def undelete_quota(self, quota_id: str) -> str:
        """
        Undelete a quota

        :type quota_id: str
        :param quota_id: Encoded quota ID.

        :rtype: str
        :return: A description of the changes, mentioning the undeleted quota.
          For example::

            "Undeleted 1 quotas: Testing-B"
        """
        url = self._make_url(quota_id, deleted=True) + "/undelete"
        return self._post(url=url, payload={"id": quota_id})
