import requests
import json


def create(database_id, name):
    entry = {
        "parent": {
            "database_id": database_id,
        },
        # our convention is to call the title of an entry "name"
        "properties": {"name": {"title": [{"text": {"content": name}}]}},
    }
    return json.dumps(entry)


class Table:
    # The implementation is based on the following two links
    # * https://developers.notion.com/docs/getting-started
    # * https://developers.notion.com/docs/working-with-databases

    def __init__(self, database_id):
        try:  # env variable fails in Jupyter notebooks
            from .._secrets import NOTION_API_KEY
        except ImportError:
            raise RuntimeError("Please run: lamin configure")

        self._base_url = "https://api.notion.com/v1/pages"
        self._headers = {
            "Authorization": f"Bearer {NOTION_API_KEY}",
            "Content-Type": "application/json",
            "Notion-Version": "2022-02-22",
        }
        self._database_id = database_id

    def __getitem__(self, id):
        response = requests.get(f"{self._base_url}/{id}", headers=self._headers)
        return response.json()

    def create(self, name):
        """Create an entry with a given name."""
        data = create(self._database_id, name)
        response = requests.post(self._base_url, headers=self._headers, data=data)
        return response.json()

    def update(self, id, **kwargs):
        """Update text fields. Other types aren't implemented right now."""
        # update the desired fields
        properties = {}
        for key, value in kwargs.items():
            properties[key] = {"rich_text": [{"text": {"content": value}}]}
        # make a json entry
        data = json.dumps({"properties": properties})
        # send patch request with updated properties
        response = requests.patch(
            f"{self._base_url}/{id}", headers=self._headers, data=data
        )
        return response.json()


class Dataset(Table):
    def __init__(self):
        super().__init__(database_id="4499c79c6e9141c8bbdcb251a24d901a")
