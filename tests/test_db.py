import lamindb as lndb


def test_create_to_load():

    lndb._admin.db.setup()

    from lamindb._configuration import user_name

    user_id = lndb._admin.db.insert_if_not_exists.user(user_name)  # type: ignore
    print(f"added user {user_id} ({user_name})")
    lndb._admin.db.insert.file("test_file.csv", source="83jf")
    for entity in lndb.schema.entities:
        print(lndb.load(entity))


if __name__ == "__main__":
    test_create_to_load()
