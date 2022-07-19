from ..Auth import Auth


class SupabaseDbClient:
    def __init__(self, auth: Auth) -> None:
        self.auth = auth

    def insert(self, table_name, params):
        data = self.auth.supabase_client.table(table_name).insert(params).execute()
        return data

    def load_rows(self, table_name, match_keys):
        statement = self.auth.supabase_client.table(table_name).select("*")
        for key, value in match_keys.items():
            statement = statement.eq(key, value)
        data = statement.execute()
        return data

    def delete_rows(self, table_name, match_keys):
        statement = self.auth.supabase_client.table(table_name).delete()
        for key, value in match_keys.items():
            statement = statement.eq(key, value)
        data = statement.execute()
        return data
