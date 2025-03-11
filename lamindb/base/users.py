user_id_cache = {}


def current_user_id() -> int:
    import lamindb_setup as ln_setup
    from lamindb_setup import settings
    from lamindb_setup._init_instance import register_user

    from lamindb.models import User

    def query_user_id():
        if ln_setup.core.django.IS_MIGRATING:
            return 1
        else:
            try:
                user_id = User.objects.get(uid=settings.user.uid).id
            except User.DoesNotExist:
                register_user(settings.user)
                user_id = User.objects.get(uid=settings.user.uid).id
            return user_id

    if settings._instance_exists:
        if settings.instance.slug not in user_id_cache:
            user_id_cache[settings.instance.slug] = query_user_id()
        return user_id_cache[settings.instance.slug]
    else:
        return query_user_id()
