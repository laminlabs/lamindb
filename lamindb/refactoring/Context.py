from lamindb.refactoring.instance import load_instance
from lamindb.refactoring.settings import ContextSettingsStore, create_settings_model


class Context:
    @staticmethod
    def get_current_user():
        context_settings = Context.__create_or_get_context_settings()
        assert context_settings["current_user_email"] is not None
        assert context_settings["current_user_secret"] is not None
        return {
            "current_user_email": context_settings["current_user_email"],
            "current_user_secret": context_settings["current_user_secret"],
        }

    @staticmethod
    def set_current_user(email: str, secret: str):
        context_settings = Context.__create_or_get_context_settings()
        context_settings["current_user_email"] = email
        context_settings["current_user_secret"] = secret

    @staticmethod
    def get_current_instance():
        context_settings = Context.__create_or_get_context_settings()
        assert context_settings["current_instance_name"] is not None
        assert context_settings["current_user_email"] is not None
        assert context_settings["current_user_secret"] is not None
        current_instance = load_instance(
            context_settings["current_instance_name"],
            context_settings["current_user_email"],
        )
        return current_instance

    @staticmethod
    def set_current_instance(instance_name: str):
        context_settings = Context.__create_or_get_context_settings()
        context_settings["current_instance_name"] = instance_name

    @staticmethod
    def __create_or_get_context_settings():
        ContextSettings = create_settings_model(ContextSettingsStore)
        settings = ContextSettings("context")
        settings_store = ContextSettingsStore()
        if not settings.exists:
            settings.setup(settings_store)
        return settings
