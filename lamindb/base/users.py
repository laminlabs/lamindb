user_id_cache = {}


def _has_write_access() -> bool:
    from django.db import connection

    with connection.cursor() as cursor:
        cursor.execute("""
            EXISTS (
                SELECT 1 FROM check_access() chk
                WHERE chk.role in ('write', 'admin')
            )
        """)
        return cursor.fetchone()


def current_user_id() -> int:
    import lamindb_setup as ln_setup
    from lamindb_setup import settings
    from lamindb_setup._init_instance import register_user

    from lamindb.errors import NoWriteAccess
    from lamindb.models import User

    def query_user_id():
        if ln_setup.core.django.IS_MIGRATING:
            return 1
        else:
            user = settings.user
            user_uid = user.uid
            try:
                user_id = User.objects.get(uid=user_uid).id
            except User.DoesNotExist:
                register_user(user)
                try:
                    user_id = User.objects.get(uid=user_uid).id
                except User.DoesNotExist as e:
                    isettings = settings.instance
                    if isettings.is_read_only_connection:
                        raise NoWriteAccess(
                            "Unable to register a new user in the instance database "
                            "because you have a read-only connection."
                        ) from e
                    if isettings._db_permissions == "jwt" and not _has_write_access():
                        raise NoWriteAccess(
                            "Unable to register a new user in the instance database "
                            "because you don't have write access to any space or registry."
                        ) from e
                    raise e
            return user_id

    if settings._instance_exists:
        slug = settings.instance.slug
        if slug not in user_id_cache:
            user_id_cache[slug] = query_user_id()
        return user_id_cache[slug]
    else:
        return query_user_id()
