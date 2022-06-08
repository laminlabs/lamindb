import lamindb as lndb


def test_create_to_load():

    lndb.db.meta.create()

    from lamindb._configuration import user_name

    user_id = lndb.db.insert_if_not_exists.user(user_name)  # type: ignore
    print(f"added user {user_id} ({user_name})")
    lndb.db.insert.file("test_file.csv", "83jf")
    for entity in lndb.db.entities:
        print(lndb.db.load(entity))


if __name__ == "__main__":
    test_create_to_load()
