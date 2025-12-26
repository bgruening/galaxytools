"""
Contains possible interactions with the Galaxy visualization
"""

from typing import (
    Any,
    TYPE_CHECKING,
)

from bioblend.galaxy.client import Client

if TYPE_CHECKING:
    from bioblend.galaxy import GalaxyInstance


class VisualClient(Client):
    module = "visualizations"

    def __init__(self, galaxy_instance: "GalaxyInstance") -> None:
        super().__init__(galaxy_instance)

    def get_visualizations(self) -> list[dict[str, Any]]:
        """
        Get the list of all visualizations.

        :rtype: list
        :return: A list of dicts with details on individual visualizations.
          For example::

            [{'dbkey': 'eschColi_K12',
              'id': 'df1c7c96fc427c2d',
              'title': 'AVTest1',
              'type': 'trackster',
              'url': '/api/visualizations/df1c7c96fc427c2d'},
             {'dbkey': 'mm9',
              'id': 'a669f50f8bf55b02',
              'title': 'Bam to Bigwig',
              'type': 'trackster',
              'url': '/api/visualizations/a669f50f8bf55b02'}]
        """
        return self._get()

    def show_visualization(self, visual_id: str) -> dict[str, Any]:
        """
        Get details of a given visualization.

        :type visual_id: str
        :param visual_id: Encoded visualization ID

        :rtype: dict
        :return: A description of the given visualization.
          For example::

            {'annotation': None,
             'dbkey': 'mm9',
             'id': '18df9134ea75e49c',
             'latest_revision': {  ... },
             'model_class': 'Visualization',
             'revisions': ['aa90649bb3ec7dcb', '20622bc6249c0c71'],
             'slug': 'visualization-for-grant-1',
             'title': 'Visualization For Grant',
             'type': 'trackster',
             'url': '/u/azaron/v/visualization-for-grant-1',
             'user_id': '21e4aed91386ca8b'}
        """
        return self._get(id=visual_id)
