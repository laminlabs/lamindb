from supabase import create_client

from lamindb.refactoring.utils.get_default_supabase_credentials import (
    get_default_supabase_credentials,
)

from .Context import Context


class Auth:
    def __init__(self, url: str = None, key: str = None) -> None:
        if not (url and key):
            url, key = get_default_supabase_credentials()
        self.supabase_client = create_client(url, key)
        self.session = None

    def sign_up(self, email: str, secret: str) -> None:
        user = self.supabase_client.auth.sign_up(email=email, password=secret)
        if not user.identities:
            raise RuntimeError(
                "It seems you already have an account. Please use login function"
            )
        if self.__is_already_waiting_for_confirmation(user):
            RuntimeWarning("Email confirmation already sended")
        return user

    def log_in(self, email: str, secret: str) -> None:
        self.log_in_supabase(email, secret)
        Context.set_current_user(email, secret)
        self.log_out_supabase()

    def log_out(self):
        Context.set_current_user(None, None)

    def log_in_supabase_current_user(self) -> None:
        current_user = Context.get_current_user()
        self.log_in_supabase(
            current_user["current_user_email"], current_user["current_user_password"]
        )

    def log_in_supabase(self, email: str, secret: str) -> None:
        try:
            self.session = self.supabase_client.auth.sign_in(
                email=email, password=secret
            )
            if self.session:
                self.supabase_client.postgrest.auth(self.session.access_token)
        except Exception:
            raise Exception("Wrong email or secret.")

    def log_out_supabase(self) -> None:
        self.supabase_client.auth.sign_out()

    def __is_already_waiting_for_confirmation(self, user) -> bool:
        diff = user.confirmation_sent_at - user.identities[0].last_sign_in_at
        return diff.total_seconds() > 0.1
